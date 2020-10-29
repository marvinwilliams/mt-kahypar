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
#include "mt-kahypar/partition/metrics.h"
#include "mt-kahypar/utils/memory_tree.h"
#include "mt-kahypar/utils/timer.h"

namespace mt_kahypar {

void BasicGreedyRefiner::initializeImpl(PartitionedHypergraph &phg) {
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

  // don't forget to set the new imbalance and km1 values in the metrics object.
  // you can ignore the cut value
  if (!_is_initialized) {
    throw std::runtime_error("Call initialize before calling refine");
  }
  LOG << "You called greedy refinement";

  Gain overall_improvement = 0;
  size_t consecutive_rounds_with_too_little_improvement = 0;
  sharedData.release_nodes = context.refinement.fm.release_nodes;
  tbb::task_group tg;
  vec<HypernodeWeight> initialPartWeights(size_t(sharedData.numParts));
  HighResClockTimepoint fm_start = std::chrono::high_resolution_clock::now();
  utils::Timer &timer = utils::Timer::instance();

  for (size_t round = 0; round < context.refinement.fm.multitry_rounds;
       ++round) { // global multi try rounds
    for (PartitionID i = 0; i < sharedData.numParts; ++i) {
      initialPartWeights[i] = phg.partWeight(i);
    }

    timer.start_timer("collect_border_nodes", "Collect Border Nodes");
    roundInitialization(phg);
    timer.stop_timer("collect_border_nodes");

    size_t num_border_nodes = sharedData.refinementNodes.unsafe_size();
    if (num_border_nodes == 0) {
      break;
    }
    size_t num_seeds = context.refinement.fm.num_seed_nodes;
    if (context.type == kahypar::ContextType::main &&
        !refinement_nodes.empty() /* n-level */
        && num_border_nodes < 20 * context.shared_memory.num_threads) {
      num_seeds = num_border_nodes / (4 * context.shared_memory.num_threads);
      num_seeds = std::min(num_seeds, context.refinement.fm.num_seed_nodes);
      num_seeds = std::max(num_seeds, 1UL);
    }

    timer.start_timer("find_moves", "Find Moves");
    sharedData.finishedTasks.store(0, std::memory_order_relaxed);
    auto task = [&](const size_t task_id) {
      auto &greedy = ets_bgf.local();
      while (sharedData.finishedTasks.load(std::memory_order_relaxed) <
                 sharedData.finishedTasksLimit &&
             greedy.findMoves(phg, task_id, num_seeds)) { /* keep running*/
      }
      sharedData.finishedTasks.fetch_add(1, std::memory_order_relaxed);
    };
    size_t num_tasks =
        std::min(num_border_nodes, context.shared_memory.num_threads);
    ASSERT(static_cast<int>(num_tasks) <=
           TBBNumaArena::instance().total_number_of_threads());
    for (size_t i = 0; i < num_tasks; ++i) {
      tg.run(std::bind(task, i));
    }
    tg.wait();
    timer.stop_timer("find_moves");

    timer.start_timer("rollback", "Rollback to Best Solution");
    HyperedgeWeight improvement = globalRollback.revertToBestPrefix<
        GainCacheStrategy::maintain_gain_cache_between_rounds>(
        phg, sharedData, initialPartWeights);
    timer.stop_timer("rollback");

    const double roundImprovementFraction =
        improvementFraction(improvement, metrics.km1 - overall_improvement);
    overall_improvement += improvement;
    if (roundImprovementFraction < context.refinement.fm.min_improvement) {
      consecutive_rounds_with_too_little_improvement++;
    } else {
      consecutive_rounds_with_too_little_improvement = 0;
    }

    HighResClockTimepoint fm_timestamp =
        std::chrono::high_resolution_clock::now();
    const double elapsed_time =
        std::chrono::duration<double>(fm_timestamp - fm_start).count();
    if (debug && context.type == kahypar::ContextType::main) {
      FMStats stats;
      for (auto &fm : ets_bgf) {
        fm.stats.merge(stats);
      }
      LOG << V(round) << V(improvement) << V(metrics::km1(phg))
          << V(metrics::imbalance(phg, context)) << V(num_border_nodes)
          << V(roundImprovementFraction) << V(elapsed_time);
//          << V(current_time_limit) << stats.serialize();
    }

    if (improvement <= 0 ||
        consecutive_rounds_with_too_little_improvement >= 2) {
      break;
    }
  }

  if (context.partition.show_memory_consumption &&
      context.partition.verbose_output &&
      context.type == kahypar::ContextType::main &&
      phg.initialNumNodes() ==
          sharedData.moveTracker.moveOrder.size() /* top level */) {
    printMemoryConsumption();
  }

  metrics.km1 -= overall_improvement;
  metrics.imbalance = metrics::imbalance(phg, context);
  ASSERT(metrics.km1 == metrics::km1(phg), V(metrics.km1)
                                               << V(metrics::km1(phg)));
  return overall_improvement > 0;

  return false;
}

/* TODO: task queue strategy via param <28-10-20, @noahares> */
void BasicGreedyRefiner::roundInitialization(PartitionedHypergraph &phg) {
  // clear border nodes
  sharedData.refinementNodes.clear();

  // iterate over all nodes and insert border nodes into task queue
  tbb::parallel_for(
      tbb::blocked_range<HypernodeID>(0, phg.initialNumNodes()),
      [&](const tbb::blocked_range<HypernodeID> &r) {
        const int task_id = tbb::this_task_arena::current_thread_index();
        ASSERT(task_id >= 0 &&
               task_id < TBBNumaArena::instance().total_number_of_threads());
        for (HypernodeID u = r.begin(); u < r.end(); ++u) {
          if (phg.nodeIsEnabled(u) && phg.isBorderNode(u)) {
            sharedData.refinementNodes.safe_push(u, task_id);
          }
        }
      });

  // shuffle task queue if requested
  if (context.refinement.fm.shuffle) {
    sharedData.refinementNodes.shuffle();
  }

  // requesting new searches activates all nodes by raising the deactivated node
  // marker also clears the array tracking search IDs in case of overflow
  sharedData.nodeTracker.requestNewSearches(
      static_cast<SearchID>(sharedData.refinementNodes.unsafe_size()));
}

} // namespace mt_kahypar
