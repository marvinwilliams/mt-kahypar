//
// Created by mlaupichler on 19.04.21.
//

#pragma once

#include <queue>
#include <functional>

/* #define L1_CACHE_LINESIZE 64 */
/* #define PAGESIZE 4096 */

#include "utils/priority_queue_factory.hpp"

//#include <tbb/concurrent_queue.h>
/* #include "mt-kahypar/definitions.h" */
#include "mt-kahypar/partition/context.h"
#include "uncontraction_group_tree.h"
#include "mt-kahypar/datastructures/async/async_common.h"
#include "depth_priority_queue.h"
#include "node_region_comparator.h"
#include "mt-kahypar/parallel/parallel_prefix_sum.h"


namespace mt_kahypar::ds
{
    using DoParallelForAllGroupsFunction = std::function<void (const ContractionGroup&)>;
    using DoParallelForAllGroupIDsFunction = std::function<void (const ContractionGroupID&)>;

    template<typename GroupHierarchy = Mandatory,
        typename RegionComparator = Mandatory>
    class ConcurrentQueueGroupPool {

    public:

        using pq_type =
          typename util::PriorityQueueFactory<unsigned long, unsigned long>::type;

        /**
         * Constructs new sequential contraction group pool. Passes ownership of the given group hierarchy to this pool.
         */
        explicit ConcurrentQueueGroupPool(std::unique_ptr<GroupHierarchy> hierarchy, const Context &context)
                : _hierarchy(hierarchy.release()),
                  _pools_by_task(context.shared_memory.num_threads),
                  _node_region_comparator(nullptr),
                  _region_similarity_retries(context.uncoarsening.node_region_similarity_retries)
                  , _picked_ids_ets(std::vector<ContractionGroupID>(context.uncoarsening.node_region_similarity_retries, invalidGroupID)),
                  _pick_candidates_ets(std::vector<pick_candidate>(context.uncoarsening.node_region_similarity_retries)),
                  _total_calls_to_pick(0),
                  _calls_to_pick_that_reached_max_retries(0),
                  _calls_to_pick_with_empty_pq(0),
                  _num_accepted_uncontractions(0),
                  _calculate_full_similarities(context.uncoarsening.region_comparison_with_full_similarities),
                  _group_pool_type(context.uncoarsening.group_pool_type),
                  _similarity_to_others_threshold(context.uncoarsening.node_region_similarity_threshold)
                {
            ASSERT(_hierarchy);
	    ASSERT(_group_pool_type == GroupPoolType::multiqueue);
	    auto pq_params = util::PriorityQueueParameters{};
#ifdef MQ_C
	    pq_params.c = MQ_C;
#endif
#ifdef MQ_STICKINESS
	    pq_params.stickiness = MQ_STICKINESS;
#endif
	    _pq = std::unique_ptr<pq_type>{new pq_type(util::create_pq<pq_type>(1'000'000, context.shared_memory.num_threads, pq_params))};
	    for (auto root : _hierarchy->roots()) {
		    _pq->push({_hierarchy->depth(root), root});
	    }

        }

        ConcurrentQueueGroupPool(const ConcurrentQueueGroupPool& other) = delete;
        ConcurrentQueueGroupPool& operator=(const ConcurrentQueueGroupPool& other) = delete;
        ConcurrentQueueGroupPool(ConcurrentQueueGroupPool&& other) = delete;
        ConcurrentQueueGroupPool& operator=(ConcurrentQueueGroupPool&& other) = delete;

        ~ConcurrentQueueGroupPool() = default;

        const GroupHierarchy* getPtrToHierarchyForQueries() {
          ASSERT(_hierarchy);
          return _hierarchy.get();
        }

        void setNodeRegionComparator(RegionComparator* comparator) {
          ASSERT(comparator);
          _node_region_comparator = comparator;
        }

        // ===== Hierarchy forwards =====

        uint32_t getNumTotal() {
          ASSERT(_hierarchy);
            return _hierarchy->getNumGroups();
        }

        HypernodeID getTotalNumUncontractions() const {
          ASSERT(_hierarchy);
          return _hierarchy->getNumContainedContractions();
        }

        size_t getVersion() const {
          ASSERT(_hierarchy);
            return _hierarchy->getVersion();
        }

        const ContractionGroup &group(const ContractionGroupID id) const {
          ASSERT(_hierarchy);
            return _hierarchy->group(id);
        }

        const HypernodeID& depth(const ContractionGroupID id) const {
          ASSERT(_hierarchy);
            return _hierarchy->depth(id);
        }

        const HypernodeID getNumberOfUncontractionsInDepth(const HypernodeID depth) const {
          ASSERT(_hierarchy);
          return _hierarchy->getNumberOfUncontractionsInDepth(depth);
        }

        bool isLastContractionGroupOfRepresentative(const ContractionGroupID groupID) const {
          ASSERT(_hierarchy);
          return _hierarchy->isLastContractionGroupOfRepresentative(groupID);
        }

        // ! Returns whether a node is initially stable for this version, i.e. whether it is not the representative in
        // ! any uncontraction in the underlying hierarchy
        bool isInitiallyStableNode(const HypernodeID hn) const {
          ASSERT(_hierarchy);
          return _hierarchy->isInitiallyStableNode(hn);
        }

        ContractionGroupIDIteratorRange successors(const ContractionGroupID id) {
          ASSERT(_hierarchy);
          return _hierarchy->successors(id);
        }

        ContractionGroupID numSuccessors(const ContractionGroupID id) {
          ASSERT(_hierarchy);
          return _hierarchy->numSuccessors(id);
        }

        // ===== Activation/Deactivation of Hypernodes =====

        bool tryToPickActiveID(ContractionGroupID& destination, const size_t task_id, typename pq_type::Handle& handle) {
          if (_calculate_full_similarities) {
            return tryToPickActiveIDWithFullSimilarities(destination, task_id, handle);
          } else {
            return tryToPickActiveIDWithEarlyBreak(destination, task_id, handle);
          }
        }

        bool tryToPickActiveIDWithFullSimilarities(ContractionGroupID& destination, const size_t task_id, typename pq_type::Handle& handle, const bool accept_based_on_similarity_to_others = true) {
          _total_calls_to_pick.add_fetch(1, std::memory_order_relaxed);

            destination = invalidGroupID;
            if (!_node_region_comparator || (_region_similarity_retries == 0)) {
              return tryPopActive(destination, handle);
            } else {
              const size_t num_edges_active_in_this_and_other_tasks = _node_region_comparator->numberOfEdgesActiveInThisTaskAndAnyOtherTask(task_id);
              std::vector<pick_candidate>& candidates = _pick_candidates_ets.local();
              ContractionGroupID best_idx = 0;
              for (uint32_t i = 0; i < _region_similarity_retries; ++i) {
                pick_candidate next_candidate;
                bool picked = tryToGetNewCandidate(next_candidate, task_id, handle, num_edges_active_in_this_and_other_tasks);
                if (!picked) {
                  // If none found in PQ and none have been previously, return false (no elements available)
                  if (i == 0) {
                    _calls_to_pick_with_empty_pq.fetch_add(1, std::memory_order_relaxed);
                    return false;
                  }
                  // If no more found in PQ but previously some ID was found, reinsert all but best and return best
                  destination = candidates[best_idx].group_id;
                  for (uint32_t j = 0; j < i; ++j) {
                    if (j != best_idx) {
                      insertActive(candidates[j].group_id, handle);
                    }
                  }
                  if (destination == invalidGroupID) {
                    ERROR("Destination invalid after picking best from so far picked when PQ popping failed!" << V(destination));
                  }
                  return true;
                }
                ASSERT(next_candidate.group_id != invalidGroupID);
                if (next_candidate.group_id == invalidGroupID) {
                  ERROR("No correct group has been picked!" << V(next_candidate.group_id));
                }
                if (isGoodCandidate(next_candidate, accept_based_on_similarity_to_others)) {
                  // If good candidate found, return it right away and reinsert others that were picked
                  destination = next_candidate.group_id;
                  for (uint32_t j = 0; j < i; ++j) {
                      insertActive(candidates[j].group_id, handle);
                  }
                  return true;
                } else {
                  // If candidate is not good, insert it into candidates and possibly update best
                  candidates[i] = next_candidate;
                  if (isBetterThanOtherCandidate(next_candidate, candidates[best_idx])) {
                    best_idx = i;
                  }
                }
              }
              // If tried as often as possible, reinsert all but best and return best
              _calls_to_pick_that_reached_max_retries.fetch_add(1, std::memory_order_relaxed);
              destination = candidates[best_idx].group_id;
              for (uint32_t j = 0; j < _region_similarity_retries; ++j) {
                if (j != best_idx) {
                  insertActive(candidates[j].group_id, handle);
                }
              }
              ASSERT(destination != invalidGroupID);
              return true;
            }
        }

        bool tryToPickActiveIDWithEarlyBreak(ContractionGroupID& destination, const size_t task_id, typename pq_type::Handle& handle) {
          _total_calls_to_pick.add_fetch(1, std::memory_order_relaxed);

          destination = invalidGroupID;
          if (!_node_region_comparator || (_region_similarity_retries == 0)) {
            return tryPopActive(destination, handle);
          } else {
            std::vector<ContractionGroupID>& pickedIDs = _picked_ids_ets.local();
            size_t best_idx = 0;
            double best_part_of_edges_before_fail = 0.0;
            for (uint32_t i = 0; i < _region_similarity_retries; ++i) {
              ContractionGroupID pickedID = invalidGroupID;
              bool picked = tryPopActive(pickedID, handle);
              if (!picked) {
                if (i == 0) {
                  _calls_to_pick_with_empty_pq.fetch_add(1, std::memory_order_relaxed);
                  return false;
                } else {
                  // If no more found in PQ and no good one was found just return the one that failed at the latest (heuristic) and reinsert the rest
                  destination = pickedIDs[best_idx];
                  for (uint32_t j = 0; j < i; ++j) {
                    if (j != best_idx) {
                      insertActive(pickedIDs[j], handle);
                    }
                  }
                }
                ASSERT(destination != invalidGroupID);
                return true;
              }
              ASSERT(pickedID != invalidGroupID);
              HypernodeID repr = group(pickedID).getRepresentative();
              double part_of_edges_seen_before_failure;
              if (_node_region_comparator->regionIsNotTooSimilarToActiveNodesWithEarlyBreak(repr, task_id, part_of_edges_seen_before_failure)) {
                // If good candidate found return it and reinsert all previously found
                destination = pickedID;
                for (uint32_t j = 0; j < i; ++j) {
                  insertActive(pickedIDs[j], handle);
                }
                return true;
              } else {
                pickedIDs[i] = pickedID;
                if (part_of_edges_seen_before_failure > best_part_of_edges_before_fail) {
                  best_idx = i;
                  best_part_of_edges_before_fail = part_of_edges_seen_before_failure;
                }
              }
            }
            // If tried as often as possible, reinsert all but the one that failed at the latest (heuristic)
            _calls_to_pick_that_reached_max_retries.fetch_add(1, std::memory_order_relaxed);
            destination = pickedIDs[best_idx];
            for (uint32_t j = 0; j < _region_similarity_retries; ++j) {
              if (j != best_idx) {
                insertActive(pickedIDs[j], handle);
              }
            }
            ASSERT(destination != invalidGroupID);
            return true;
          }
        }

        void activate(const ContractionGroupID id, typename pq_type::Handle& handle) {
            insertActive(id, handle);
        }

        // ! Mark a given id as accepted, i.e. mark that a thread will commence working on uncoarsening the group and
        // ! that it will not be reinserted into the pool again.
        void markAccepted(const ContractionGroupID id) {
          ASSERT(_hierarchy);
//          if (!_use_multiqueue) {
//            ASSERT(_dpq_active_ids);
//            _dpq_active_ids->increment_finished(_hierarchy->depth(id));
//          }
          _num_accepted_uncontractions.fetch_add(group(id).size(), std::memory_order_relaxed);
        }

        size_t getNumAccepted() const {
          return _num_accepted_uncontractions.load(std::memory_order_relaxed);
        }

        // ===== Parallel Iteration Convenience Methods =====

        void doParallelForAllGroups(const DoParallelForAllGroupsFunction& f) const {
            tbb::parallel_for(all(),[&](BlockedGroupIDIterator& range) {
                for (ContractionGroupID id = range.begin(); id != range.end(); ++id) {
                    f(group(id));
                }
            });
        }

        void doParallelForAllGroupIDs(const DoParallelForAllGroupIDsFunction& f) const {
            tbb::parallel_for(all(),[&](BlockedGroupIDIterator& range) {
                for (ContractionGroupID id = range.begin(); id != range.end(); ++id) {
                    f(id);
                }
            });
        }

        size_t getTotalCallsToPick() const {
          return _total_calls_to_pick.load(std::memory_order_relaxed);
        }

        size_t getCallsToPickThatReachedMaxRetries() const {
          return _calls_to_pick_that_reached_max_retries.load(std::memory_order_relaxed);
        }

        size_t getCallsToPickWithEmptyPQ() const {
          return _calls_to_pick_with_empty_pq.load(std::memory_order_relaxed);
        }

        void memoryConsumption(utils::MemoryTreeNode* parent) const {
          ASSERT(parent);

          utils::MemoryTreeNode* hierarchy_node = parent->addChild("Uncontraction Hierarchy");
          _hierarchy->memoryConsumption(hierarchy_node);
//          if (!_use_multiqueue) {
//            ASSERT(_dpq_active_ids);
//            utils::MemoryTreeNode* depth_pq_node = parent->addChild("Depth Priority Queue");
//            _dpq_active_ids->memoryConsumption(depth_pq_node);
//          }
          if (_group_pool_type == GroupPoolType::thread_local_pools) {
            size_t total_size_of_pools = 0;
            for (size_t i = 0; i < _pools_by_task.size(); ++i) {
              total_size_of_pools += _pools_by_task[i]->size_in_bytes();
            }
            utils::MemoryTreeNode* pools_per_task_node = parent->addChild("Active Groups Pools");
            pools_per_task_node->updateSize(total_size_of_pools);
          }
        }

        bool taskFinished() const {
            return allAccepted();
        }

        bool allAccepted() const {
          size_t totalNumUncontractions = getTotalNumUncontractions();
          size_t numAcceptedUncontractions = _num_accepted_uncontractions.load(std::memory_order_relaxed);
          return numAcceptedUncontractions == totalNumUncontractions;
        }

	pq_type::Handle getPQHandle() {
		return _pq->get_handle();
	}

    private:

        bool poolEmptyForTask(const size_t task_id) const {
          ASSERT(_group_pool_type == GroupPoolType::thread_local_pools);
          ASSERT(task_id < _pools_by_task.size());
          return _pools_by_task[task_id]->empty();
        }

        BlockedGroupIDIterator all() const {
          ASSERT(_hierarchy);
          return _hierarchy->all();
        }

//        bool tryInsertActive(ContractionGroupID id) {
//          ASSERT(_hierarchy);
////          _active_ids.push(id, _hierarchy->depth(id));
////          return true;
//          bool locked = _queue_lock.tryLock();
//          if (!locked) return false;
//          _active_ids.push(id);
//          _queue_lock.unlock();
//          return true;
//        }


        bool tryPopActive(ContractionGroupID& destination, typename pq_type::Handle& handle) {
            typename pq_type::value_type ret;
            bool success = handle.try_extract_top(ret);
            if (success) {
              destination = ret.second;
            }
            return success;
        }

        void insertActive(ContractionGroupID id, typename pq_type::Handle& handle) {
            ASSERT(_hierarchy);
            handle.push({_hierarchy->depth(id), id});
        }

        struct pick_candidate {
            ContractionGroupID group_id = invalidGroupID;
            double similarity_to_calling_task = 0.0;
            double similarity_to_other_tasks = 1.0;
        };

        bool tryToGetNewCandidate(pick_candidate& candidate, const size_t task_id, typename pq_type::Handle& handle, const HyperedgeID num_active_in_this_and_other_tasks) {
          bool picked = tryPopActive(candidate.group_id, handle);
          if (!picked) return false;

          HypernodeID repr = group(candidate.group_id).getRepresentative();
          auto [sim_to_this_task, sim_to_other_tasks] = _node_region_comparator->regionSimilarityToThisTaskAndOtherTasks(repr, task_id, num_active_in_this_and_other_tasks);
          candidate.similarity_to_calling_task = sim_to_this_task;
          candidate.similarity_to_other_tasks = sim_to_other_tasks;
          return true;
        }

        bool isGoodCandidate(const pick_candidate& candidate, const bool accept_based_only_on_similarity_to_others) const {
          return (candidate.similarity_to_other_tasks < _similarity_to_others_threshold + CMP_EPSILON
            && (accept_based_only_on_similarity_to_others || candidate.similarity_to_calling_task > _similarity_to_calling_task_threshold - CMP_EPSILON));
        }

        // ! Returns true iff this_candidate is better than other_candidate
        bool isBetterThanOtherCandidate(const pick_candidate& this_candidate, const pick_candidate& other_candidate) const {

          // If both are below threshold for other tasks, the one with greater similarity to this task is better
          if (this_candidate.similarity_to_other_tasks < _similarity_to_others_threshold + CMP_EPSILON
              && other_candidate.similarity_to_other_tasks < _similarity_to_others_threshold + CMP_EPSILON) {
            return this_candidate.similarity_to_calling_task > other_candidate.similarity_to_calling_task;
          }

          // Otherwise the one with the lower similarity to other tasks is better
          return (this_candidate.similarity_to_other_tasks < other_candidate.similarity_to_other_tasks
            || (this_candidate.similarity_to_other_tasks <= other_candidate.similarity_to_other_tasks + CMP_EPSILON
              && this_candidate.similarity_to_calling_task > other_candidate.similarity_to_calling_task));
        }

        std::unique_ptr<GroupHierarchy> _hierarchy;

	std::unique_ptr<pq_type> _pq;
        // todo: work stealing for thread local pools
        std::vector<std::unique_ptr<NoDownsizeIntegralTypeVector<ContractionGroupID>>> _pools_by_task;

        RegionComparator* _node_region_comparator;
        const size_t _region_similarity_retries;
        tbb::enumerable_thread_specific<std::vector<ContractionGroupID>> _picked_ids_ets;
        tbb::enumerable_thread_specific<std::vector<pick_candidate>> _pick_candidates_ets;

        CAtomic<size_t> _total_calls_to_pick;
        CAtomic<size_t> _calls_to_pick_that_reached_max_retries;
        CAtomic<size_t> _calls_to_pick_with_empty_pq;

        CAtomic<size_t> _num_accepted_uncontractions;

        const bool _calculate_full_similarities;
        const GroupPoolType _group_pool_type;

        const double _similarity_to_others_threshold;
        const double _similarity_to_calling_task_threshold = 0.05;

        static constexpr double CMP_EPSILON = 1.0e-30;


    };

//    using TreeGroupPool = ConcurrentQueueGroupPool<UncontractionGroupTree, Hypergraph>;
//    using VersionedPoolVector = parallel::scalable_vector<std::unique_ptr<TreeGroupPool>>;


} // namespace mt_kahypar::ds
