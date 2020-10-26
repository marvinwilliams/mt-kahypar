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

#include "mt-kahypar/partition/refinement/greedy/greedy_refiner.h"
#include "mt-kahypar/datastructures/hypergraph_common.h"

namespace mt_kahypar {

void BasicGreedyRefiner::initializeImpl(PartitionedHypergraph &phg) {
  /* TODO: timer needed? <17-10-20, @noahares> */
  if (!phg.isGainCacheInitialized()) {
    phg.initializeGainCache();
  }
  _is_initialized = true;
}

/**
 *
 * @param phg The partition and hypergraph to refine
 * @param refinement_nodes Used for n-level uncoarsening. Ignore for now. It
 * will be empty
 * @param metrics stores current km1 and imbalance
 * @return whether the refinement improved the objective function
 */
bool BasicGreedyRefiner::refineImpl(
    PartitionedHypergraph &phg,
    const parallel::scalable_vector<HypernodeID> &refinement_nodes,
    kahypar::Metrics &metrics, double) {

  // implement the refinement

  // don't forget to set the new imbalance and km1 values in the metrics object.
  // you can ignore the cut value
  if (!_is_initialized) {
    throw std::runtime_error("Call initialize on greedy before calling refine");
  }

  Gain total_gain = 0;
  vec<HyperedgeWeight> initial_part_weights(size_t(shared_data.num_parts));

  // setup partition weights for balance constraints
  for (PartitionID i = 0; i < shared_data.num_parts; ++i) {
    initial_part_weights[i] = phg.partWeight(i);
  }

  initBorderVertices(phg, refinement_nodes);

  /* TODO: keep refinement nodes up to date <24-10-20, @noahares> */
  /* TODO: what happens if highest gain move is not allowed due to balance? what
   * if this holds for all refinement nodes? <24-10-20, @noahares> */
  for (HypernodeID v : refinement_nodes) {
    auto[best_to, gain] = findBestToPartition(phg, v);
    shared_data.target_part[v] = best_to;
    shared_data.gain_buckets.insert(gain, v);
    shared_data.best_gain = std::max(gain, shared_data.best_gain);
  }

  auto move_candidates =
      shared_data.gain_buckets.getBucket(shared_data.best_gain);

  auto node_to_move =
      std::find_if(move_candidates.begin(), move_candidates.end(),
                   [&](const HypernodeID move_candidate) {
                     phg.changeNodePartWithGainCacheUpdate(
                         move_candidate, phg.partID(move_candidate),
                         shared_data.target_part[move_candidate]);
                   });

  if (node_to_move == move_candidates.end()) {
    LOG << "No move with maximum gain possible";
    /* TODO: what to do in this case? move on to next best gain (where to store
     * it?) <25-10-20, @noahares> */
  } else {
    updateNeighbors(phg, *node_to_move);
  }

  LOG << "You called greedy refinement";

  return false;
}

void BasicGreedyRefiner::initBorderVertices(
    PartitionedHypergraph &phg,
    const parallel::scalable_vector<HypernodeID> &refinement_nodes) {
  shared_data.refinement_nodes.clear();
  if (refinement_nodes.empty()) {
    tbb::parallel_for(
        tbb::blocked_range<HypernodeID>(0, phg.initialNumNodes()),
        [&](const tbb::blocked_range<HypernodeID> &r) {
          const int task_id = tbb::this_task_arena::current_thread_index();
          ASSERT(task_id >= 0 &&
                 task_id < TBBNumaArena::instance().total_number_of_threads());
          for (HypernodeID u = r.begin(); u < r.end(); ++u) {
            if (phg.nodeIsEnabled(u) && phg.isBorderNode(u)) {
              shared_data.refinement_nodes.safe_push(u, task_id);
            }
          }
        });
  }
}

std::pair<PartitionID, Gain>
BasicGreedyRefiner::findBestToPartition(const PartitionedHypergraph &phg,
                                        const HypernodeID u) {
  const HypernodeWeight u_weight = phg.nodeWeight(u);
  const PartitionID from = phg.partID(u);
  const HypernodeWeight from_weight = phg.partWeight(from);
  PartitionID to = kInvalidPartition;
  HypernodeWeight to_penalty = std::numeric_limits<HyperedgeWeight>::max();
  HypernodeWeight to_best_weight = from_weight - u_weight;
  const auto improvement = [&](HypernodeWeight to_weight,
                               HypernodeWeight penalty, PartitionID i) {
    return ((penalty < to_penalty ||
             (penalty == to_penalty && to_weight < to_best_weight)) &&
            to_weight + u_weight <= _context.partition.max_part_weights[i]);
  };
  for (PartitionID i = 0; i < phg.k(); ++i) {
    if (i == from) {
      continue;
    }
    const HypernodeWeight to_weight = phg.partWeight(i);
    const HypernodeWeight penalty = phg.moveToPenalty(u, i);
    if (improvement(to_weight, penalty, i)) {
      to_penalty = penalty;
      to = i;
      to_best_weight = to_weight;
    }
  }

  const Gain gain = to != kInvalidPartition
                        ? phg.moveFromBenefit(u) - to_penalty
                        : std::numeric_limits<HyperedgeWeight>::min();
  return std::make_pair(to, gain);
}

void BasicGreedyRefiner::updateNeighbors(PartitionedHypergraph &phg,
                                         const HypernodeID u) {
  for (auto e : phg.incidentEdges(u)) {
    for (auto v : phg.pins(e)) {
      /* TODO: update gains (and part weights) <25-10-20, @noahares> */
    }
  }
}

} // namespace mt_kahypar
