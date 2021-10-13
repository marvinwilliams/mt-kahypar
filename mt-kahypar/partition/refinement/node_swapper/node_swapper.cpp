/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2019 Tobias Heuer <tobias.heuer@kit.edu>
 * Copyright (C) 2020 Lars Gottesbüren <lars.gottesbueren@kit.edu>
 *
 * KaHyPar is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * KaHyPar is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with KaHyPar.  If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/

#include "mt-kahypar/partition/refinement/node_swapper/node_swapper.h"

#include <queue>

#include "tbb/concurrent_vector.h"
#include "tbb/parallel_sort.h"

#include "mt-kahypar/parallel/atomic_wrapper.h"
#include "mt-kahypar/utils/randomize.h"


namespace mt_kahypar {

#define NOOP_FUNC [] (const HyperedgeID, const HyperedgeWeight, const HypernodeID, const HypernodeID, const HypernodeID) { }


namespace {

class VertexMoveComparator {
public:
  bool operator()(const NodeSwapper::VertexMove& lhs, const NodeSwapper::VertexMove& rhs) {
    return lhs.gain < rhs.gain;
  }
};

// Simple concurrent priority queue.
// Our concurrent priority queue contains p (number of threads) priority queues.
// push and pop selects a random priority queue and protect writes to it via a spin lock.
class ConcurrentPQ {

using VertexMove = typename NodeSwapper::VertexMove;
using PQ = std::priority_queue<VertexMove, std::vector<VertexMove>, VertexMoveComparator>;

public:
  ConcurrentPQ(const size_t num_threads) :
    _num_threads(num_threads),
    _pq_lock(num_threads),
    _pqs(num_threads) { }

  void push(VertexMove&& move) {
    const size_t i = selectPQ();
    _pq_lock[i].lock();
    _pqs[i].emplace(std::move(move));
    _pq_lock[i].unlock();
  }

  VertexMove pop() {
    const size_t i = selectPQ();
    _pq_lock[i].lock();
    if ( !_pqs[i].empty() ) {
      VertexMove move = _pqs[i].top();
      _pqs[i].pop();
      _pq_lock[i].unlock();
      return move;
    }
    _pq_lock[i].unlock();
    return VertexMove { kInvalidHypernode, kInvalidPartition, kInvalidPartition, kInvalidGain };
  }

private:
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE size_t selectPQ() {
    return utils::Randomize::instance().getRandomInt(
      0, static_cast<int>(_num_threads) - 1, sched_getcpu());
  }

  const size_t _num_threads;
  vec<SpinLock> _pq_lock;
  vec<PQ> _pqs;
};

NodeSwapper::VertexMove computeBestTargetBlock(const PartitionedHypergraph& phg,
                                               const HypernodeID hn,
                                               const PartitionID k,
                                               const bool force_move) {
  const PartitionID from = phg.partID(hn);
  NodeSwapper::VertexMove move { hn, from, from, force_move ? std::numeric_limits<Gain>::min() : 0 };
  for ( PartitionID to = 0; to < k; ++to ) {
    if ( from != to ) {
      const Gain gain = phg.km1Gain(hn, from, to);
      if ( gain > move.gain || ( ( !force_move || gain > 0 ) && gain == move.gain &&
          utils::Randomize::instance().flipCoin(sched_getcpu()) )) {
        move.gain = gain;
        move.to = to;
      }
    }
  }
  return move;
}

} // namespace

HyperedgeWeight NodeSwapper::refine() {
  if ( !_hg.isGainCacheInitialized() ) {
    _hg.initializeGainCache();
  }
  ASSERT(_hg.checkTrackedPartitionInformation());

  const HyperedgeWeight initial_km1 = metrics::km1(_hg);
  HyperedgeWeight delta = optimisitcHighDegreeNodeMoving(initial_km1);

  // We perform several rounds until the overall improvement is below 0.25%
  HyperedgeWeight last_round_km1 = initial_km1;
  double last_round_reduction = 2.0;
  while ( last_round_reduction > 0.0025 ) {
    delta += nodeSwappingRound();
    const HyperedgeWeight current_km1 = initial_km1 - delta;
    last_round_reduction = static_cast<double>(last_round_km1) / current_km1 - 1.0;
    last_round_km1 = current_km1;
  }

  return -delta;
}

HyperedgeWeight NodeSwapper::optimisitcHighDegreeNodeMoving(const HyperedgeWeight current_km1) {
  CAtomic<HyperedgeWeight> delta(0);
  _in_queue.reset();
  for ( vec<VertexMove>& moves : _local_moves ) {
    moves.clear();
  }

  tbb::enumerable_thread_specific<vec<VertexMove>> local_moves;
  vec<tbb::enumerable_thread_specific<HyperedgeWeight>> local_block_volumes(
    _context.partition.k, tbb::enumerable_thread_specific<HyperedgeWeight>(0));
  _hg.doParallelForAllEdges([&](const HyperedgeID& he) {
    const HyperedgeWeight edge_weight = _hg.edgeWeight(he);
    for ( const PartitionID block : _hg.connectivitySet(he) ) {
      local_block_volumes[block].local() += edge_weight;
    }
  });
  vec<HyperedgeWeight> block_volumes(_context.partition.k, 0);
  PartitionID target = kInvalidPartition;
  HyperedgeWeight heaviest_weight = 0;
  for ( PartitionID i = 0; i < _context.partition.k; ++i ) {
    block_volumes[i] = local_block_volumes[i].combine(std::plus<HyperedgeWeight>());
    if ( block_volumes[i] > heaviest_weight ) {
      target = i;
      heaviest_weight = block_volumes[i];
    }
  }

  _hg.doParallelForAllNodes([&](const HypernodeID& hn) {
    const PartitionID from = _hg.partID(hn);
    HyperedgeWeight cut_contribution = 0;
    if ( from != target ) {
      for ( const HyperedgeID he : _hg.incidentEdges(hn) ) {
        if ( _hg.pinCountInPart(he, target) > 0 ) {
          cut_contribution += _hg.edgeWeight(he);
        }
      }
    }

    if ( static_cast<double>(cut_contribution) / current_km1 > 0.01 ) {
      HyperedgeWeight gain = 0;
      auto delta_gain_func = [&](const HyperedgeID he,
                                 const HyperedgeWeight edge_weight,
                                 const HypernodeID edge_size,
                                 const HypernodeID pin_count_in_from_part_after,
                                 const HypernodeID pin_count_in_to_part_after) {
        gain -= km1Delta(he, edge_weight, edge_size,
          pin_count_in_from_part_after, pin_count_in_to_part_after);
      };

      _hg.changeNodePartWithGainCacheUpdate(hn, from, target,
        std::numeric_limits<HypernodeWeight>::max(), [] { }, delta_gain_func);
      _local_moves.local().emplace_back(VertexMove { hn, from, target, gain });

      _in_queue.set(hn);
      delta += gain;
    }
  });

  if ( _hg.partWeight(target) > _context.partition.max_part_weights[target] ) {
    ConcurrentPQ pq(_context.shared_memory.num_threads);
    _hg.doParallelForAllNodes([&](const HypernodeID& hn) {
      const PartitionID from = _hg.partID(hn);
      if ( from == target && !_in_queue[hn] ) {
        pq.push(computeBestTargetBlock(_hg, hn, _context.partition.k, true));
      }
    });

    CAtomic<HyperedgeWeight> target_weight(_hg.partWeight(target));
    tbb::parallel_for(0UL, _context.shared_memory.num_threads, [&](const size_t) {
      while ( true ) {
        VertexMove move = pq.pop();
        if ( move.hn != kInvalidHypernode && target_weight > _context.partition.max_part_weights[target] ) {

          if ( move.gain != _hg.km1Gain(move.hn, move.from, move.to) ) {
            move.gain = _hg.km1Gain(move.hn, move.from, move.to);
            pq.push(std::move(move));
            continue;
          }

          HyperedgeWeight new_target_weight = target_weight.fetch_sub(_hg.nodeWeight(move.hn));
          if ( new_target_weight > _context.partition.max_part_weights[target] ) {
            HyperedgeWeight gain = 0;
            auto delta_gain_func = [&](const HyperedgeID he,
                                      const HyperedgeWeight edge_weight,
                                      const HypernodeID edge_size,
                                      const HypernodeID pin_count_in_from_part_after,
                                      const HypernodeID pin_count_in_to_part_after) {
              gain -= km1Delta(he, edge_weight, edge_size,
                pin_count_in_from_part_after, pin_count_in_to_part_after);
            };

            if ( _hg.changeNodePartWithGainCacheUpdate(move.hn, move.from, move.to,
              _context.partition.max_part_weights[move.to], [] { }, delta_gain_func) ) {
              _local_moves.local().emplace_back(std::move(move));
              _in_queue.set(move.hn);
              delta += gain;
            } else {
              target_weight += _hg.nodeWeight(move.hn);
            }
          }
        } else {
          break;
        }
      }
    });
  }

  if ( _hg.partWeight(target) < _context.partition.max_part_weights[target] ) {
    ConcurrentPQ pq(_context.shared_memory.num_threads);
    _hg.doParallelForAllNodes([&](const HypernodeID& hn) {
      const PartitionID from = _hg.partID(hn);
      if ( from != target && !_in_queue[hn] ) {
        VertexMove move = computeBestTargetBlock(_hg, hn, _context.partition.k, false);
        if ( move.to == target ) {
          pq.push(std::move(move));
        }
      }
    });

    tbb::parallel_for(0UL, _context.shared_memory.num_threads, [&](const size_t) {
      while ( true ) {
        VertexMove move = pq.pop();
        if ( move.hn != kInvalidHypernode ) {

          if ( move.gain != _hg.km1Gain(move.hn, move.from, move.to) ) {
            move.gain = _hg.km1Gain(move.hn, move.from, move.to);
            pq.push(std::move(move));
            continue;
          }

          HyperedgeWeight gain = 0;
          auto delta_gain_func = [&](const HyperedgeID he,
                                    const HyperedgeWeight edge_weight,
                                    const HypernodeID edge_size,
                                    const HypernodeID pin_count_in_from_part_after,
                                    const HypernodeID pin_count_in_to_part_after) {
            gain -= km1Delta(he, edge_weight, edge_size,
              pin_count_in_from_part_after, pin_count_in_to_part_after);
          };

          if ( _hg.changeNodePartWithGainCacheUpdate(move.hn, move.from, move.to,
            _context.partition.max_part_weights[move.to], [] { }, delta_gain_func) ) {
            _local_moves.local().emplace_back(std::move(move));
            _in_queue.set(move.hn);
            delta += gain;
          }
        } else {
          break;
        }
      }
    });
  }

  delta += nodeSwappingRound(true, false);

  if ( delta < 0 ) {
    for ( const vec<VertexMove>& moves : _local_moves ) {
      tbb::parallel_for(0UL, moves.size(), [&](const size_t i) {
        const VertexMove& move = moves[i];
        ASSERT(_hg.partID(move.hn) == move.to);
        _hg.changeNodePartWithGainCacheUpdate(move.hn, move.to, move.from);
        _hg.recomputeMoveFromBenefit(move.hn);
      });
    }
    delta.store(0, std::memory_order_relaxed);
  }

  return delta.load();
}

HyperedgeWeight NodeSwapper::nodeSwappingRound(const bool collect_moves, const bool reset_queue) {
  CAtomic<HyperedgeWeight> delta(0);
  if ( reset_queue ) {
    _in_queue.reset();
  }

  // First, we compute for each node its preferred block that improves the objective function
  // and insert it into the corresponding priority queue.
  vec<ConcurrentPQ> to_pq(_context.partition.k, ConcurrentPQ(_context.shared_memory.num_threads));
  _hg.doParallelForAllNodes([&](const HypernodeID& hn) {
    const PartitionID from = _hg.partID(hn);
    VertexMove best_move = computeBestTargetBlock(_hg, hn, _context.partition.k, false);
    if ( from != best_move.to && !_in_queue[hn] ) {
      to_pq[best_move.to].push(std::move(best_move));
      _in_queue.set(hn);
    }
  });

  // Second, for each node u not contained in one of the priority queue,
  // we try to swap it with the highest gain node that want to move into
  // the block of node u.
  utils::Randomize::instance().parallelShuffleVector(_nodes, 0UL, _nodes.size());
  tbb::parallel_for(0UL, _nodes.size(), [&](const size_t i) {
    const HypernodeID u = _nodes[i];
    if ( _hg.nodeIsEnabled(u) && !_in_queue[u] ) {
      const PartitionID u_from = _hg.partID(u);
      VertexMove move = to_pq[u_from].pop();
      while ( move.hn != kInvalidHypernode &&
              move.gain != _hg.km1Gain(move.hn, move.from, u_from) ) {
        move.gain = _hg.km1Gain(move.hn, move.from, u_from);
        if ( move.gain > 0 ) {
          to_pq[u_from].push(std::move(move));
        }
        move = to_pq[u_from].pop();
      }

      if ( move.hn != kInvalidHypernode ) {
        const HypernodeID v = move.hn;
        const PartitionID v_from = move.from;
        const Gain u_gain = _hg.km1Gain(u, u_from, v_from);

        // Note that the actual swap gain can differ, if u and v share
        // some common edges. However, we double check the estimated gain
        // with attributed gain and if the swap worsen the solution quality
        // we revert it immediatly.
        const Gain estimated_swap_gain = move.gain + u_gain;
        if ( estimated_swap_gain > 0 ) {
          Gain real_swap_gain = 0;
          auto delta_gain_func = [&](const HyperedgeID he,
                                    const HyperedgeWeight edge_weight,
                                    const HypernodeID edge_size,
                                    const HypernodeID pin_count_in_from_part_after,
                                    const HypernodeID pin_count_in_to_part_after) {
            real_swap_gain -= km1Delta(he, edge_weight, edge_size,
              pin_count_in_from_part_after, pin_count_in_to_part_after);
          };

          // Perform swap
          _hg.changeNodePartWithGainCacheUpdate(u, u_from, v_from,
            std::numeric_limits<HypernodeWeight>::max(), [] { }, delta_gain_func);
          _hg.changeNodePartWithGainCacheUpdate(v, v_from, u_from,
            std::numeric_limits<HypernodeWeight>::max(), [] { }, delta_gain_func);

          if ( real_swap_gain < 0 ) {
            // Revert swap, if it has worsen the solution quality
            _hg.changeNodePartWithGainCacheUpdate(u, v_from, u_from,
              std::numeric_limits<HypernodeWeight>::max(), [] { }, delta_gain_func);
            _hg.changeNodePartWithGainCacheUpdate(v, u_from, v_from,
              std::numeric_limits<HypernodeWeight>::max(), [] { }, delta_gain_func);
            move.gain = _hg.km1Gain(move.hn, v_from, u_from);
            to_pq[u_from].push(std::move(move));
          } else if ( collect_moves ) {
            _local_moves.local().emplace_back(VertexMove { u, u_from, v_from, 0 });
            _local_moves.local().emplace_back(VertexMove { v, v_from, u_from, 0 });
          }

          delta += real_swap_gain;
        } else {
          to_pq[u_from].push(std::move(move));
        }
      }
    }
  });

  _hg.doParallelForAllNodes([&](const HypernodeID& hn) {
    _hg.recomputeMoveFromBenefit(hn);
  });

  return delta.load();
}

} // namespace mt_kahypar