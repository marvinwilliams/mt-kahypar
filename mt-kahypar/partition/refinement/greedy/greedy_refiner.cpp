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
#include <tbb/parallel_sort.h>

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

  unused(refinement_nodes);
  Gain overall_improvement = 0;
  tbb::task_group tg;
  sharedData.release_nodes = context.refinement.greedy.release_nodes;
  HighResClockTimepoint greedy_start =
      std::chrono::high_resolution_clock::now();
  utils::Timer &timer = utils::Timer::instance();

  for (size_t round = 0; round < context.refinement.greedy.multitry_rounds;
       ++round) { // global multi try rounds

    timer.start_timer("collect_border_nodes", "Collect Border Nodes");
    roundInitialization(phg, context.refinement.greedy.assignment_strategy);
    timer.stop_timer("collect_border_nodes");

    size_t num_border_nodes = _greedy_shared_data.numRefinementNodes();
    if (num_border_nodes == 0) {
      break;
    }
    timer.start_timer("find_moves", "Find Moves");
    sharedData.finishedTasks.store(0, std::memory_order_relaxed);
    auto task = [&](const size_t task_id) {
      auto &greedy = ets_bgf.local();
      greedy.findMoves(phg, task_id);
      _greedy_shared_data.hold_barrier.lowerSize();
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

    HighResClockTimepoint greedy_timestamp =
        std::chrono::high_resolution_clock::now();
    const double elapsed_time =
        std::chrono::duration<double>(greedy_timestamp - greedy_start).count();
    FMStats stats;
    Gain improvement = 0;
    for (auto &greedy : ets_bgf) {
      greedy.stats.merge(stats);
      improvement += greedy.getGain();
      greedy.reset();
    }
    if (context.partition.verbose_output) {
      LOG << V(round) << V(improvement) << V(metrics::km1(phg))
          << V(metrics::imbalance(phg, context)) << V(num_border_nodes)
          << V(elapsed_time) << stats.serialize();
    }

    overall_improvement += improvement;
    tbb::parallel_for(MoveID(0), sharedData.moveTracker.numPerformedMoves(),
                      [&](MoveID move_id) {
                        phg.recomputeMoveFromBenefit(
                            sharedData.moveTracker.moveOrder[move_id].node);
                      });
    sharedData.moveTracker.reset();
    if (improvement <= 0)
      break;
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
}

void BasicGreedyRefiner::roundInitialization(
    PartitionedHypergraph &phg, GreedyAssignmentStrategy assignment_strategy) {
  // clear border nodes
  for (auto &v : _greedy_shared_data.refinement_nodes) {
    v.clear();
  }

  // clear message queues
  for (auto &mq : _greedy_shared_data.messages) {
    mq.clear();
  }

  _greedy_shared_data.hold_barrier.reset(context.shared_memory.num_threads);

  auto random_assignment = [&](const tbb::blocked_range<HypernodeID> &r) {
    const int task_id = tbb::this_task_arena::current_thread_index();
    ASSERT(task_id >= 0 &&
           task_id < TBBNumaArena::instance().total_number_of_threads());
    for (HypernodeID u = r.begin(); u < r.end(); ++u) {
      if (phg.nodeIsEnabled(u) && phg.isBorderNode(u)) {
        _greedy_shared_data.refinement_nodes[task_id].push_back(u);
      }
    }
  };

  switch (assignment_strategy) {
  case GreedyAssignmentStrategy::static_assignement:
    staticAssignment(phg);
    break;
  case GreedyAssignmentStrategy::random_assignment:
    tbb::parallel_for(tbb::blocked_range<HypernodeID>(0, phg.initialNumNodes()),
                      random_assignment);
    break;
  case GreedyAssignmentStrategy::partition_assignment:
    partitionAssignment(phg);
    break;
  default:
    throw std::runtime_error("no work distribution strategy provided");
  }

  // shuffle task queue if requested
  if (context.refinement.greedy.shuffle) {
    tbb::parallel_for_each(_greedy_shared_data.refinement_nodes,
                           [](vec<HypernodeID> q) {
                             utils::Randomize::instance().shuffleVector(q);
                           });
  }

  // requesting new searches activates all nodes by raising the deactivated node
  // marker also clears the array tracking search IDs in case of overflow
  sharedData.nodeTracker.requestNewSearches(
      static_cast<SearchID>(_greedy_shared_data.refinement_nodes.size()));
}

void BasicGreedyRefiner::staticAssignment(PartitionedHypergraph &phg) {

  tbb::enumerable_thread_specific<vec<HypernodeID>> ets_border_nodes;

  // thread local border node calculation
  tbb::parallel_for(tbb::blocked_range<HypernodeID>(0, phg.initialNumNodes()),
                    [&](const tbb::blocked_range<HypernodeID> &r) {
                      auto &tl_border_nodes = ets_border_nodes.local();
                      for (HypernodeID u = r.begin(); u < r.end(); ++u) {
                        if (phg.nodeIsEnabled(u) && phg.isBorderNode(u)) {
                          tl_border_nodes.push_back(u);
                        }
                      }
                    });
  // combine thread local border nodes
  vec<HypernodeID> border_nodes;
  for (const auto &tl_border_nodes : ets_border_nodes) {
    border_nodes.insert(border_nodes.end(), tl_border_nodes.begin(),
                        tl_border_nodes.end());
  }
  size_t nodes_per_thread =
      (border_nodes.size() / context.shared_memory.num_threads) + 1;
  tbb::task_group tg;
  // distribute border nodes to threads
  auto task = [&](const auto thread_id) {
    auto begin = border_nodes.begin() + thread_id * nodes_per_thread;
    auto end = std::min(begin + nodes_per_thread, border_nodes.end());
    _greedy_shared_data.refinement_nodes[thread_id].insert(
        _greedy_shared_data.refinement_nodes[thread_id].end(), begin, end);
  };
  for (size_t i = 0; i < context.shared_memory.num_threads; ++i) {
    tg.run(std::bind(task, i));
  }
  tg.wait();
}

void BasicGreedyRefiner::partitionAssignment(PartitionedHypergraph &phg) {

  tbb::enumerable_thread_specific<vec<vec<HypernodeID>>> ets_border_nodes;

  // thread local border node calculation
  tbb::parallel_for(tbb::blocked_range<HypernodeID>(0, phg.initialNumNodes()),
                    [&](const tbb::blocked_range<HypernodeID> &r) {
                      auto &tl_border_nodes = ets_border_nodes.local();
                      tl_border_nodes.resize(context.partition.k);
                      for (HypernodeID u = r.begin(); u < r.end(); ++u) {
                        if (phg.nodeIsEnabled(u) && phg.isBorderNode(u)) {
                          tl_border_nodes[phg.partID(u)].push_back(u);
                        }
                      }
                    });

  vec<std::pair<size_t, PartitionID>> part_sizes_with_id(context.partition.k);
  vec<vec<HypernodeID>> border_nodes(context.partition.k);

  // make pairs of number of border nodes in a partition and the partition id.
  // Also combine the found border nodes of each partition into one bucket
  // list
  tbb::parallel_for(PartitionID(0), context.partition.k, [&](const auto i) {
    part_sizes_with_id[i] = std::make_pair(0, i);
  });
  for (const auto &tl_border_nodes : ets_border_nodes) {
    tbb::parallel_for(PartitionID(0), context.partition.k, [&](const auto i) {
      part_sizes_with_id[i].first += tl_border_nodes[i].size();
      border_nodes[i].insert(border_nodes[i].end(), tl_border_nodes[i].begin(),
                             tl_border_nodes[i].end());
    });
  }

  // sort pairs by number of border nodes descending
  std::sort(part_sizes_with_id.begin(), part_sizes_with_id.end(),
            std::greater<>());

  size_t sum =
      std::accumulate(part_sizes_with_id.begin(), part_sizes_with_id.end(), 0,
                      [](auto a, const auto &b) { return a + b.first; });

  // optimal number of border nodes for each thread
  size_t target_size = (sum / context.shared_memory.num_threads) + 1;

  // distribute border nodes equally onto threads
  /* TODO: can bin packing be made prettier? <19-11-20, @noahares> */
  vec<vec<PartitionID>> partitions_of_thread(context.shared_memory.num_threads);
  size_t index = 0;
  for (const auto &nodes : part_sizes_with_id) {
    size_t start_index = index;
    int remaining_space = target_size -
                          _greedy_shared_data.refinement_nodes[index].size() +
                          nodes.first;
    int best_remaining_space = remaining_space;
    int best_index = index;
    bool all_bins_checked = false;
    while (remaining_space < 0 && !all_bins_checked) {
      index = (index + 1) % context.shared_memory.num_threads;
      remaining_space = target_size -
                        _greedy_shared_data.refinement_nodes[index].size() +
                        nodes.first;
      best_index = remaining_space > best_remaining_space ? index : best_index;
      best_remaining_space = std::max(best_remaining_space, remaining_space);
      if (index == start_index) {
        all_bins_checked = true;
      }
    }
    partitions_of_thread[best_index].push_back(nodes.second);
    index = (index + 1) % context.shared_memory.num_threads;
  }

  // actually copy hypernode ids to final assignment
  tbb::task_group tg;
  auto task = [&](const auto i) {
    for (auto id : partitions_of_thread[i]) {
      _greedy_shared_data.refinement_nodes[i].insert(
          _greedy_shared_data.refinement_nodes[i].end(),
          border_nodes[id].begin(), border_nodes[id].end());
    }
  };
  for (size_t i = 0; i < context.shared_memory.num_threads; ++i) {
    tg.run(std::bind(task, i));
  }
  tg.wait();

  ASSERT(sum == _greedy_shared_data.numRefinementNodes());
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
