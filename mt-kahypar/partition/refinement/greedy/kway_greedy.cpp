/*******************************************************************************
 * This file is part of MT-KaHyPar.
 *
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

#include "mt-kahypar/partition/refinement/greedy/kway_greedy.h"

namespace mt_kahypar {

void KWayGreedy::findMoves(PartitionedHypergraph &phg, size_t task_id) {
  thisSearch = ++sharedData.nodeTracker.highestActiveSearchID;
  _task_id = task_id;

  for (HypernodeID v : _greedy_shared_data.refinement_nodes[task_id]) {
    if (sharedData.nodeTracker.tryAcquireNode(v, thisSearch)) {
      fm_strategy.insertIntoPQ(phg, v, 0);
    }
  }
  fm_strategy.updatePQs(phg);

  if (runStats.pushes > 0) {
    internalFindMoves(phg);
    /* TODO: sync message queues after finishing ? <21-11-20, @noahares> */
  }
}

template <typename Partition>
MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE std::pair<PartitionID, HypernodeWeight>
heaviestPartAndWeight(const Partition &partition) {
  PartitionID p = kInvalidPartition;
  HypernodeWeight w = std::numeric_limits<HypernodeWeight>::min();
  for (PartitionID i = 0; i < partition.k(); ++i) {
    if (partition.partWeight(i) > w) {
      w = partition.partWeight(i);
      p = i;
    }
  }
  return std::make_pair(p, w);
}

template <typename PHG>
MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE void
KWayGreedy::updateNeighbors(PHG &phg, const Move &move) {
  // Note: In theory we should acquire/update all neighbors. It just turned out
  // that this works fine Actually: only vertices incident to edges with gain
  // changes can become new boundary vertices. Vertices that already were
  // boundary vertices, can still be considered later since they are in the task
  // queue
  // --> actually not that bad
  for (const auto &egu : edgesWithGainChanges) {
    HyperedgeID e = egu.e;
    if (phg.edgeSize(e) < context.partition.ignore_hyperedge_size_threshold) {
      for (HypernodeID v : phg.pins(e)) {
        if (neighborDeduplicator[v] != deduplicationTime) {
          SearchID searchOfV = sharedData.nodeTracker.searchOfNode[v].load(
              std::memory_order_acq_rel);
          if (searchOfV == thisSearch) {
            fm_strategy.updateGain(phg, v, move);
          } else if (context.refinement.greedy.sync_with_mq && searchOfV != 0 &&
                     searchOfV > sharedData.nodeTracker.deactivatedNodeMarker) {
            // send hypernode id to responsible threads message queue
            SearchID v_index =
                searchOfV - sharedData.nodeTracker.deactivatedNodeMarker - 1;
            SearchID this_index =
                thisSearch - sharedData.nodeTracker.deactivatedNodeMarker - 1;
            SearchID num_threads = context.shared_memory.num_threads;
            ASSERT(v_index * num_threads + this_index <
                   static_cast<SearchID>(_greedy_shared_data.messages.size()));
            _greedy_shared_data.messages[v_index * num_threads + this_index]
                .push_back(v);
          }
          neighborDeduplicator[v] = deduplicationTime;
        }
      }
    }
  }
  edgesWithGainChanges.clear();

  updateNeighborDeduplicator();
}

void KWayGreedy::internalFindMoves(PartitionedHypergraph &phg) {
  Move move;

  auto delta_func = [&](const HyperedgeID he, const HyperedgeWeight edge_weight,
                        const HypernodeID,
                        const HypernodeID pin_count_in_from_part_after,
                        const HypernodeID pin_count_in_to_part_after) {
    // Gains of the pins of a hyperedge can only change in the following
    // situations.
    if (pin_count_in_from_part_after == 0 ||
        pin_count_in_from_part_after == 1 || pin_count_in_to_part_after == 1 ||
        pin_count_in_to_part_after == 2) {
      edgesWithGainChanges.push_back(
          {he, pin_count_in_from_part_after, pin_count_in_to_part_after});
    }

    _gain += (pin_count_in_to_part_after == 1 ? -edge_weight : 0) +
             (pin_count_in_from_part_after == 0 ? edge_weight : 0);
  };
  auto delta_func_revert =
      [&](const HyperedgeID, const HyperedgeWeight edge_weight,
          const HypernodeID, const HypernodeID pin_count_in_from_part_after,
          const HypernodeID pin_count_in_to_part_after) {
        _gain += (pin_count_in_to_part_after == 1 ? -edge_weight : 0) +
                 (pin_count_in_from_part_after == 0 ? edge_weight : 0);
      };

  size_t bestImprovementIndex = 0;
  Gain estimatedImprovement = 0;
  Gain bestImprovement = 0;
  Gain lastImprovement = std::numeric_limits<Gain>::max();

  HypernodeWeight heaviestPartWeight = 0;
  HypernodeWeight fromWeight = 0, toWeight = 0;

  while (lastImprovement > 0 &&
         sharedData.finishedTasks.load(std::memory_order_relaxed) <
             sharedData.finishedTasksLimit) {

    if (_local_moves_since_sync >=
        context.refinement.greedy.num_moves_before_sync) {
      if (context.refinement.greedy.sync_with_mq) {
        syncMessageQueues(phg);
      } else {
        syncAllLocalNodes(phg);
      }
    }

    if (!fm_strategy.findNextMoveNoRetry(phg, move))
      break;

    sharedData.nodeTracker.deactivateNode(move.node, thisSearch);
    MoveID move_id = std::numeric_limits<MoveID>::max();
    bool moved = false;
    Gain delta_before = _gain;
    if (move.to != kInvalidPartition) {
      heaviestPartWeight = heaviestPartAndWeight(phg).second;
      fromWeight = phg.partWeight(move.from);
      toWeight = phg.partWeight(move.to);
      moved = phg.changeNodePart(
          move.node, move.from, move.to,
          context.partition.max_part_weights[move.to],
          [&] { move_id = sharedData.moveTracker.insertMove(move); },
          delta_func);
    }

    if (moved) {
      const Gain move_delta = _gain - delta_before;
      const bool accept_move =
          (move_delta == move.gain && move_delta >= 0) || move_delta >= 0;
//      LOG << "diff in deltas" << V(move_delta) << V(move.gain);
      if (accept_move) {
        for (const auto &egu : edgesWithGainChanges) {
          // perform directly on phg and not abstract through fm_strategy
          // because we dont have a delta_phg
          phg.gainCacheUpdate(egu.e, phg.edgeWeight(egu.e), move.from,
                              egu.pin_count_in_from_part_after, move.to,
                              egu.pin_count_in_to_part_after);
        }
        runStats.moves++;
        lastImprovement = move_delta;
        estimatedImprovement += move_delta;
        localMoves.emplace_back(move, move_id);
        _local_moves_since_sync++;
        const bool improved_km1 = move_delta > 0;
        const bool improved_balance_less_equal_km1 =
            fromWeight == heaviestPartWeight &&
            toWeight + phg.nodeWeight(move.node) < heaviestPartWeight;

        if (improved_km1 || improved_balance_less_equal_km1) {
          bestImprovement = estimatedImprovement;
          bestImprovementIndex = localMoves.size();
        }
        updateNeighbors(phg, move);
      } else {
        // implicit max_to_weight = std::numeric_limits<>::max()
        phg.changeNodePart(move.node, move.to, move.from, delta_func_revert);
        sharedData.moveTracker.invalidateMove(move_id);
        edgesWithGainChanges.clear();
      }
    }

    fm_strategy.updatePQs(phg);
  }

  runStats.estimated_improvement = bestImprovement;
  fm_strategy.clearPQs(bestImprovementIndex);
  runStats.merge(stats);
}

void KWayGreedy::syncMessageQueues(PartitionedHypergraph &phg) {
  _greedy_shared_data.hold_barrier.aquire();
  SearchID this_index =
      thisSearch - sharedData.nodeTracker.deactivatedNodeMarker - 1;
  size_t num_threads = context.shared_memory.num_threads;
  auto mq_begin =
      _greedy_shared_data.messages.begin() + this_index * num_threads;
  auto mq_end = mq_begin + num_threads;
  ASSERT(mq_end <= _greedy_shared_data.messages.end());
  ASSERT((mq_begin + this_index)->empty());
  std::for_each(mq_begin, mq_end, [&](auto &mq) {
    for (const auto v : mq) {
      // use deduplicator to prevent uneeded pq updates
      if (neighborDeduplicator[v] != deduplicationTime &&
          !sharedData.nodeTracker.isLocked(v)) {
        fm_strategy.updateGainFromOtherSearch(phg, v);
        neighborDeduplicator[v] = deduplicationTime;
      }
    }
    fm_strategy.updatePQs(phg);
    mq.clear();
  });
  updateNeighborDeduplicator();
  _local_moves_since_sync = 0;
  _greedy_shared_data.hold_barrier.release();
}

void KWayGreedy::syncAllLocalNodes(PartitionedHypergraph &phg) {

  for (const auto v : _greedy_shared_data.refinement_nodes[_task_id]) {
    if (!sharedData.nodeTracker.isLocked(v)) {
      fm_strategy.updateGainFromOtherSearch(phg, v);
    }
  }
  fm_strategy.updatePQs(phg);
}

void KWayGreedy::memoryConsumption(utils::MemoryTreeNode *parent) const {
  ASSERT(parent);

  utils::MemoryTreeNode *localized_fm_node = parent->addChild("k-Way Greedy");

  utils::MemoryTreeNode *deduplicator_node =
      localized_fm_node->addChild("Deduplicator");
  deduplicator_node->updateSize(neighborDeduplicator.capacity() *
                                sizeof(HypernodeID));
  utils::MemoryTreeNode *edges_to_activate_node =
      localized_fm_node->addChild("edgesWithGainChanges");
  edges_to_activate_node->updateSize(edgesWithGainChanges.capacity() *
                                     sizeof(EdgeGainUpdate));

  utils::MemoryTreeNode *local_moves_node = parent->addChild("Local FM Moves");
  local_moves_node->updateSize(localMoves.capacity() *
                               sizeof(std::pair<Move, MoveID>));

  fm_strategy.memoryConsumption(localized_fm_node);
  // TODO fm_strategy.memoryConsumptiom(..)
  /*
  utils::MemoryTreeNode* block_pq_node = localized_fm_node->addChild("Block
  PQ"); block_pq_node->updateSize(blockPQ.size_in_bytes());
  utils::MemoryTreeNode* vertex_pq_node = localized_fm_node->addChild("Vertex
  PQ"); for ( const VertexPriorityQueue& pq : vertexPQs ) {
    vertex_pq_node->updateSize(pq.size_in_bytes());
  }
   */
}

} // namespace mt_kahypar
