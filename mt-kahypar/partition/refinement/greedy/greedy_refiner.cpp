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
  sharedData.release_nodes = context.refinement.greedy.release_nodes;
  HighResClockTimepoint greedy_start =
      std::chrono::high_resolution_clock::now();
  utils::Timer &timer = utils::Timer::instance();

  for (size_t round = 0; round < context.refinement.greedy.multitry_rounds;
       ++round) { // global multi try rounds

    // clear message queues
    for (auto i : _messages) {
      for (auto queue : i) {
        parallel::scalable_queue<HyperedgeID>().swap(queue);
      }
    }

    timer.start_timer("collect_border_nodes", "Collect Border Nodes");
    determineRefinementNodes(phg);
    timer.stop_timer("collect_border_nodes");

    size_t num_border_nodes = _refinement_nodes.size();
    if (num_border_nodes == 0) {
      break;
    }
    timer.start_timer("find_moves", "Find Moves");
    sharedData.finishedTasks.store(0, std::memory_order_relaxed);
    /* task groups seem to have an advantage only if tasks may be appended
     * later on, but here this is not the case, so parallel_for makes a good
     * choice */
    /* TODO: back to task queue <07-11-20, @noahares> */
    tbb::parallel_for(
        tbb::blocked_range<size_t>(0, _refinement_nodes.size()),
        [&](const tbb::blocked_range<size_t> &r) {
          auto &greedy = ets_bgf.local();
          std::vector<HypernodeID> local_nodes(
              _refinement_nodes.begin() + r.begin(),
              _refinement_nodes.begin() + r.end());
          while (sharedData.finishedTasks.load(std::memory_order_relaxed) <
                     sharedData.finishedTasksLimit &&
                 greedy.findMoves(phg, local_nodes)) {
          }
          sharedData.finishedTasks.fetch_add(1, std::memory_order_relaxed);
        },
        /* TODO: check for partitioners <05-11-20, @noahares> */
        tbb::static_partitioner()); // no work stealing
    timer.stop_timer("find_moves");

    HighResClockTimepoint greedy_timestamp =
        std::chrono::high_resolution_clock::now();
    const double elapsed_time =
        std::chrono::duration<double>(greedy_timestamp - greedy_start).count();
    FMStats stats;
    Gain improvement = 0;
    for (auto &greedy : ets_bgf) {
      greedy.stats.merge(stats);
      improvement += greedy.getGain();
    }
    LOG << V(round) << V(improvement) << V(metrics::km1(phg))
        << V(metrics::imbalance(phg, context)) << V(num_border_nodes)
        << V(elapsed_time) << stats.serialize();

    if (improvement <= 0)
      break;

    overall_improvement += improvement;
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
} // namespace mt_kahypar

/* TODO: when removing the work queue round init can be done easier and static
 * <02-11-20, @noahares> */
void BasicGreedyRefiner::roundInitialization(
    PartitionedHypergraph &phg, GreedyAssigmentStrategy assignment_strategy) {
  // clear border nodes
  sharedData.refinementNodes.clear();
  CAtomic<int> task_id_static(0);

  /* TODO: other possible static assignment? since isBorderNode consumes the
   * main computation time add_fetch should not produce much overhead <01-11-20,
   * @noahares> */
  auto static_assignment = [&](const tbb::blocked_range<HypernodeID> &r) {
    for (HypernodeID u = r.begin(); u < r.end(); ++u) {
      if (phg.nodeIsEnabled(u) && phg.isBorderNode(u)) {
        const int local_task_id =
            task_id_static.add_fetch(1, std::memory_order_relaxed) %
            TBBNumaArena::instance().total_number_of_threads();
        sharedData.refinementNodes.safe_push(u, local_task_id);
      }
    }
  };

  auto random_assignment = [&](const tbb::blocked_range<HypernodeID> &r) {
    const int task_id = tbb::this_task_arena::current_thread_index();
    ASSERT(task_id >= 0 &&
           task_id < TBBNumaArena::instance().total_number_of_threads());
    for (HypernodeID u = r.begin(); u < r.end(); ++u) {
      if (phg.nodeIsEnabled(u) && phg.isBorderNode(u)) {
        sharedData.refinementNodes.safe_push(u, task_id);
      }
    }
  };

  auto partition_assignment = [&](const tbb::blocked_range<HypernodeID> &r) {
    for (HypernodeID u = r.begin(); u < r.end(); ++u) {
      if (phg.nodeIsEnabled(u) && phg.isBorderNode(u)) {
        const PartitionID task_id =
            phg.partID(u) % TBBNumaArena::instance().total_number_of_threads();
        ASSERT(task_id >= 0 &&
               task_id < TBBNumaArena::instance().total_number_of_threads());
        /* TODO: mt-methis better for load balancing <2020-11-02, @noahares> */
        sharedData.refinementNodes.safe_push(u, task_id);
      }
    }
  };

  switch (assignment_strategy) {
  case GreedyAssigmentStrategy::static_assignement:
    LOG << "Round initialized with static assignment strategy";
    tbb::parallel_for(tbb::blocked_range<HypernodeID>(0, phg.initialNumNodes()),
                      static_assignment);
    break;
  case GreedyAssigmentStrategy::random_assignment:
    LOG << "Round initialized with random assignment strategy";
    tbb::parallel_for(tbb::blocked_range<HypernodeID>(0, phg.initialNumNodes()),
                      random_assignment);
    break;
  case GreedyAssigmentStrategy::partition_assignment:
    LOG << "Round initialized with partition assignment strategy";
    tbb::parallel_for(tbb::blocked_range<HypernodeID>(0, phg.initialNumNodes()),
                      partition_assignment);
    break;
  default:
    std::runtime_error("no work distribution strategy provided");
  }

  // shuffle task queue if requested
  if (context.refinement.greedy.shuffle) {
    sharedData.refinementNodes.shuffle();
  }

  // requesting new searches activates all nodes by raising the deactivated node
  // marker also clears the array tracking search IDs in case of overflow
  sharedData.nodeTracker.requestNewSearches(
      static_cast<SearchID>(sharedData.refinementNodes.unsafe_size()));
}

void BasicGreedyRefiner::determineRefinementNodes(PartitionedHypergraph &phg) {
  _refinement_nodes.clear();
  tbb::parallel_for(tbb::blocked_range<HypernodeID>(0, phg.initialNumNodes()),
                    [&](const tbb::blocked_range<HypernodeID> &r) {
                      for (HypernodeID u = r.begin(); u < r.end(); ++u) {
                        if (phg.nodeIsEnabled(u) && phg.isBorderNode(u)) {
                          _refinement_nodes.push_back(u);
                        }
                      }
                    });
}

void BasicGreedyRefiner::printMemoryConsumption() {
  utils::MemoryTreeNode greedy_memory("k-Way Greedy",
                                      utils::OutputType::MEGABYTE);

  for (const auto &greedy : ets_bgf) {
    greedy.memoryConsumption(&greedy_memory);
  }
  sharedData.memoryConsumption(&greedy_memory);
  greedy_memory.finalize();

  //    LOG << BOLD << "\n FM Memory Consumption" << END;
  //    LOG << greedy_memory;
}

} // namespace mt_kahypar
