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

    // implement the refinement

    // don't forget to set the new imbalance and km1 values in the metrics object. you can ignore the cut value

    LOG << "You called greedy refinement";

    return false;
  }


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
