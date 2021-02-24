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

#include "mt-kahypar/partition/refinement/fm/localized_kway_fm_core.h"

namespace mt_kahypar {

  template<typename FMStrategy>
  bool LocalizedKWayFM<FMStrategy>::findMoves(PartitionedHypergraph& phg, size_t taskID, size_t numSeeds) {
    localMoves.clear();
    touched_edges.clear();
    move_edges_begin.clear();
    delayed_gain_updates.clear();
    // persistent searchIDs during one round
    if (thisSearch == 0) {
      thisSearch = ++sharedData.nodeTracker.highestActiveSearchID;
    }
    ASSERT(thisSearch - sharedData.nodeTracker.deactivatedNodeMarker <= context.shared_memory.num_threads);

    if (context.refinement.fm.random_assignment) {
      auto seeds = sharedData.shared_refinement_nodes.try_pop(numSeeds);
      if (seeds) {
        for (HypernodeID u : *seeds) {
          SearchID previousSearchOfSeedNode = sharedData.nodeTracker.searchOfNode[u].load(std::memory_order_relaxed);
          if (sharedData.nodeTracker.tryAcquireNode(u, thisSearch)) {
            fm_strategy.insertIntoPQ(phg, u, previousSearchOfSeedNode);
          }
        }
      }
    } else {
      HypernodeID seedNode;
      while (runStats.pushes < numSeeds && sharedData.refinementNodes.try_pop(seedNode, taskID)) {
        SearchID previousSearchOfSeedNode = sharedData.nodeTracker.searchOfNode[seedNode].load(std::memory_order_relaxed);
        if (sharedData.nodeTracker.tryAcquireNode(seedNode, thisSearch)) {
          fm_strategy.insertIntoPQ(phg, seedNode, previousSearchOfSeedNode);
        }
      }
    }
    fm_strategy.updatePQs(phg);

    if (runStats.pushes > 0) {
      if (!context.refinement.fm.perform_moves_global
          && deltaPhg.combinedMemoryConsumption() > sharedData.deltaMemoryLimitPerThread) {
        sharedData.deltaExceededMemoryConstraints = true;
      }

      if (sharedData.deltaExceededMemoryConstraints) {
        deltaPhg.dropMemory();
      }

      if (context.refinement.fm.perform_moves_global || sharedData.deltaExceededMemoryConstraints) {
        internalFindMoves<false>(phg);
      } else {
        deltaPhg.clear();
        deltaPhg.setPartitionedHypergraph(&phg);
        internalFindMoves<true>(phg);
      }
      if (context.refinement.fm.sync_with_mq) {
        clearMessageQueues();
      }
      return true;
    } else {
      if (context.refinement.fm.sync_with_mq) {
        clearMessageQueues();
      }
      return false;
    }
  }

  template<typename Partition>
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE std::pair<PartitionID, HypernodeWeight>
  heaviestPartAndWeight(const Partition& partition) {
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

  template<typename FMStrategy>
  template<typename PHG>
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  void LocalizedKWayFM<FMStrategy>::acquireOrUpdateNeighbors(PHG& phg, const Move& move) {
    // Note: In theory we should acquire/update all neighbors. It just turned out that this works fine
    // Actually: only vertices incident to edges with gain changes can become new boundary vertices.
    // Vertices that already were boundary vertices, can still be considered later since they are in the task queue
    // --> actually not that bad
    if (context.refinement.fm.prevent_expensive_gain_updates
        || context.refinement.fm.delay_expensive_gain_updates) {
      move_edges_begin.push_back(touched_edges.size());
    }
    for (HyperedgeID e : edgesWithGainChanges) {
      if ((context.refinement.fm.prevent_expensive_gain_updates
           || context.refinement.fm.delay_expensive_gain_updates) &&
          phg.edgeSize(e) >= context.refinement.fm.large_he_threshold) {
          touched_edges.push_back(e);
      }
      if (phg.edgeSize(e) < context.partition.ignore_hyperedge_size_threshold) {
        for (HypernodeID v : phg.pins(e)) {
          if (neighborDeduplicator[v] != deduplicationTime) {
            SearchID searchOfV = sharedData.nodeTracker.searchOfNode[v].load(std::memory_order_acq_rel);
            if (searchOfV == thisSearch) {
              fm_strategy.updateGain(phg, v, move);
            } else if (sharedData.nodeTracker.tryAcquireNode(v, thisSearch)) {
              fm_strategy.insertIntoPQ(phg, v, searchOfV);
            } else if (context.refinement.fm.sync_with_mq &&
                searchOfV > sharedData.nodeTracker.deactivatedNodeMarker) {
              // send hypernode id to responsible threads message queue
              SearchID v_index = searchOfV - sharedData.nodeTracker.deactivatedNodeMarker - 1;
              SearchID this_index = thisSearch - sharedData.nodeTracker.deactivatedNodeMarker - 1;
              SearchID num_threads = context.shared_memory.num_threads;
              ASSERT(v_index * num_threads + this_index <
                  static_cast<SearchID>(sharedData.messages.size()));
              while(!sharedData.messages[v_index * num_threads + this_index]
                .try_write(v)) {};
            }
            neighborDeduplicator[v] = deduplicationTime;
          }
        }
      }
    }
    edgesWithGainChanges.clear();

    updateNeighborDeduplicator();
  }


  template<typename FMStrategy>
  template<bool use_delta>
  void LocalizedKWayFM<FMStrategy>::internalFindMoves(PartitionedHypergraph& phg) {
    StopRule stopRule(phg.initialNumNodes());
    Move move;

    auto delta_func = [&](const HyperedgeID he,
                          const HyperedgeWeight edge_weight,
                          const HypernodeID,
                          const HypernodeID pin_count_in_from_part_after,
                          const HypernodeID pin_count_in_to_part_after) {
      // Gains of the pins of a hyperedge can only change in the following situations.
      if (pin_count_in_from_part_after == 0 || pin_count_in_from_part_after == 1 ||
          pin_count_in_to_part_after == 1 || pin_count_in_to_part_after == 2) {
        if (context.refinement.fm.delay_expensive_gain_updates && moveForbidden(phg, move)) {
          delayed_gain_updates.push_back({move.node, he, pin_count_in_from_part_after, pin_count_in_to_part_after});
          return;
        }
        edgesWithGainChanges.push_back(he);
      }

      if constexpr (use_delta) {
        fm_strategy.deltaGainUpdates(deltaPhg, he, edge_weight, move.from, pin_count_in_from_part_after,
                                     move.to, pin_count_in_to_part_after);
        // prevent gain updates here
      } else {
        fm_strategy.deltaGainUpdates(phg, he, edge_weight, move.from, pin_count_in_from_part_after,
                                     move.to, pin_count_in_to_part_after);
      }

    };

    // we can almost make this function take a generic partitioned hypergraph
    // we would have to add the success func to the interface of DeltaPhg (and then ignore it there...)
    // and do the local rollback outside this function


    size_t bestImprovementIndex = 0;
    Gain estimatedImprovement = 0;
    Gain bestImprovement = 0;

    HypernodeWeight heaviestPartWeight = 0;
    HypernodeWeight fromWeight = 0, toWeight = 0;

    while (!stopRule.searchShouldStop()
           && sharedData.finishedTasks.load(std::memory_order_relaxed) < sharedData.finishedTasksLimit) {

      if (context.refinement.fm.sync_with_mq && _local_moves_since_sync >=
          context.refinement.fm.num_moves_before_sync) {
        if constexpr (use_delta) {
          syncMessageQueues(deltaPhg);
        } else {
          syncMessageQueues(phg);
        }
      }

      if constexpr (use_delta) {
        if (!fm_strategy.findNextMove(deltaPhg, move)) break;
      } else {
        if (!fm_strategy.findNextMove(phg, move)) break;
      }

      if (context.refinement.fm.prevent_expensive_gain_updates && moveForbidden(phg, move)) {
        /* TODO: really deactive node?
           theoretically the node should still be able to move.
           The problem is that we would need to reinsert it into the PQ,
           but most likely it will get the same (forbidden) target block with a potential high gain
           and therefore will get extracted when searching for the next move <12-02-21, @noahares> */
        sharedData.nodeTracker.deactivateNode(move.node, thisSearch);
        if constexpr (use_delta) {
          fm_strategy.updatePQs(deltaPhg);
        } else {
          fm_strategy.updatePQs(phg);
        }
        continue;
      }

      sharedData.nodeTracker.deactivateNode(move.node, thisSearch);
      MoveID move_id = std::numeric_limits<MoveID>::max();
      bool moved = false;
      if (move.to != kInvalidPartition) {
        /*ASSERT([&]{ return !moveForbidden(phg, move); }());*/
        if constexpr (use_delta) {
          heaviestPartWeight = heaviestPartAndWeight(deltaPhg).second;
          fromWeight = deltaPhg.partWeight(move.from);
          toWeight = deltaPhg.partWeight(move.to);
          moved = deltaPhg.changeNodePart(move.node, move.from, move.to,
                                          context.partition.max_part_weights[move.to], delta_func);
        } else {
          heaviestPartWeight = heaviestPartAndWeight(phg).second;
          fromWeight = phg.partWeight(move.from);
          toWeight = phg.partWeight(move.to);
          moved = phg.changeNodePart(move.node, move.from, move.to,
                                     context.partition.max_part_weights[move.to],
                                     [&] { move_id = sharedData.moveTracker.insertMove(move); }, delta_func);
        }
      }

      if (moved) {
        runStats.moves++;
        estimatedImprovement += move.gain;
        localMoves.emplace_back(move, move_id);
        if (context.refinement.fm.sync_with_mq) {
          _local_moves_since_sync++;
        }
        stopRule.update(move.gain);
        const bool improved_km1 = estimatedImprovement > bestImprovement;
        const bool improved_balance_less_equal_km1 = estimatedImprovement >= bestImprovement
                                                     && fromWeight == heaviestPartWeight
                                                     && toWeight + phg.nodeWeight(move.node) < heaviestPartWeight;

        if (improved_km1 || improved_balance_less_equal_km1) {
          stopRule.reset();
          bestImprovement = estimatedImprovement;
          bestImprovementIndex = localMoves.size();
        }

        if constexpr (use_delta) {
          acquireOrUpdateNeighbors(deltaPhg, move);
        } else {
          acquireOrUpdateNeighbors(phg, move);
        }
      }

      if constexpr (use_delta) {
        fm_strategy.updatePQs(deltaPhg);
      } else {
        fm_strategy.updatePQs(phg);
      }

    }

    if constexpr (use_delta) {
      std::tie(bestImprovement, bestImprovementIndex) =
              applyMovesOnGlobalHypergraph(phg, bestImprovementIndex, bestImprovement);
    } else {
      revertToBestLocalPrefix(phg, bestImprovementIndex);
    }

    if ((context.refinement.fm.prevent_expensive_gain_updates || context.refinement.fm.delay_expensive_gain_updates)
        && !touched_edges.empty()) {
      updateExpensiveMoveRevertCounter(bestImprovementIndex);
    }

    runStats.estimated_improvement = bestImprovement;
    fm_strategy.clearPQs(bestImprovementIndex);
    runStats.merge(stats);
  }


  template<typename FMStrategy>
  std::pair<Gain, size_t> LocalizedKWayFM<FMStrategy>::applyMovesOnGlobalHypergraph(
          PartitionedHypergraph& phg,
          const size_t bestGainIndex,
          const Gain bestEstimatedImprovement) {
    // TODO find better variable names!

    Gain estimatedImprovement = 0;
    Gain lastGain = 0;

    auto delta_gain_func = [&](const HyperedgeID he,
                               const HyperedgeWeight edge_weight,
                               const HypernodeID edge_size,
                               const HypernodeID pin_count_in_from_part_after,
                               const HypernodeID pin_count_in_to_part_after) {
      lastGain += km1Delta(he, edge_weight, edge_size,
                           pin_count_in_from_part_after, pin_count_in_to_part_after);
    };

    // Apply move sequence to original hypergraph and update gain values
    Gain bestImprovement = 0;
    size_t bestIndex = 0;
    for (size_t i = 0; i < bestGainIndex; ++i) {
      Move& local_move = localMoves[i].first;
      MoveID& move_id = localMoves[i].second;
      lastGain = 0;

      if constexpr (FMStrategy::uses_gain_cache) {
        phg.changeNodePartWithGainCacheUpdate(local_move.node, local_move.from, local_move.to,
                                              std::numeric_limits<HypernodeWeight>::max(),
                                              [&] { move_id = sharedData.moveTracker.insertMove(local_move); },
                                              delta_gain_func);
      } else {
        phg.changeNodePart(local_move.node, local_move.from, local_move.to,
                           std::numeric_limits<HypernodeWeight>::max(),
                           [&] { move_id = sharedData.moveTracker.insertMove(local_move); },
                           delta_gain_func);
      }

      lastGain = -lastGain; // delta func yields negative sum of improvements, i.e. negative values mean improvements
      estimatedImprovement += lastGain;
      ASSERT(move_id != std::numeric_limits<MoveID>::max());
      Move& global_move = sharedData.moveTracker.getMove(move_id);
      global_move.gain = lastGain; // Update gain value based on hypergraph delta
      if (estimatedImprovement >= bestImprovement) {  // TODO also incorporate balance into this?
        bestImprovement = estimatedImprovement;
        bestIndex = i;
      }
    }

    runStats.local_reverts += localMoves.size() - bestGainIndex;
    if (bestIndex != bestGainIndex) {
      runStats.best_prefix_mismatch++;
    }

    // Kind of double rollback, if gain values are not correct
    if (estimatedImprovement < 0) {
      // always using the if-branch gave similar results
      runStats.local_reverts += bestGainIndex - bestIndex + 1;
      for (size_t i = bestIndex + 1; i < bestGainIndex; ++i) {
        Move& m = sharedData.moveTracker.getMove(localMoves[i].second);

        if constexpr (FMStrategy::uses_gain_cache) {
          phg.changeNodePartWithGainCacheUpdate(m.node, m.to, m.from);
        } else {
          phg.changeNodePart(m.node, m.to, m.from);
        }

        sharedData.moveTracker.invalidateMove(m);
      }
      return std::make_pair(bestImprovement, bestIndex);
    } else {
      return std::make_pair(bestEstimatedImprovement, bestGainIndex);
    }
  }

  template<typename FMStrategy>
  void LocalizedKWayFM<FMStrategy>::revertToBestLocalPrefix(PartitionedHypergraph& phg, size_t bestGainIndex) {
    runStats.local_reverts += localMoves.size() - bestGainIndex;
    size_t next_revert = localMoves.size();
    while (next_revert > bestGainIndex) {
      Move& m = sharedData.moveTracker.getMove(localMoves[next_revert - 1].second);
      if constexpr (FMStrategy::uses_gain_cache) {
        phg.changeNodePartWithGainCacheUpdate(m.node, m.to, m.from);
      } else {
        phg.changeNodePart(m.node, m.to, m.from);
      }
      sharedData.moveTracker.invalidateMove(m);
      --next_revert;
    }
    if (context.refinement.fm.delay_expensive_gain_updates && !delayed_gain_updates.empty()) {
      if (bestGainIndex == 0) {
        return;
      }
      auto first_move_to_keep = localMoves.begin() + bestGainIndex - 1;
      auto next_gain_update_to_apply = delayed_gain_updates.begin();
      for (auto i = localMoves.begin(); i <= first_move_to_keep
           && next_gain_update_to_apply != delayed_gain_updates.end(); ++i) {
        Move &m = i->first;
        while (next_gain_update_to_apply != delayed_gain_updates.end()
               && next_gain_update_to_apply->node == m.node) {
          auto &d = next_gain_update_to_apply;
          phg.gainCacheUpdate(d->edge, phg.edgeWeight(d->edge), m.from, d->pin_count_in_from_part_after, m.to, d->pin_count_in_to_part_after);
          ++next_gain_update_to_apply;
        }
      }
      delayed_gain_updates.clear();
    }
  }

  template<typename FMStrategy>
  void LocalizedKWayFM<FMStrategy>::updateExpensiveMoveRevertCounter(size_t bestGainIndex) {
    move_edges_begin.push_back(touched_edges.size());
    auto next_expensive_move = move_edges_begin.begin() + bestGainIndex;
    auto next_move_to_revert = localMoves.begin() + bestGainIndex;
    for (auto i = next_move_to_revert; i < localMoves.end(); ++i) {
      Move& m = i->first;
      auto begin = touched_edges.begin() + *next_expensive_move;
      auto end = touched_edges.begin() + *(++next_expensive_move);
      ASSERT(begin <= touched_edges.end() && end <= touched_edges.end());
      for (auto e = begin; e < end; ++e) {
        ASSERT(*e < sharedData.num_edges_up_to.size());
        size_t index = sharedData.num_edges_up_to[*e];
        ASSERT(index < sharedData.num_large_he);
        ASSERT(index * context.partition.k + m.to < sharedData.forbidden_move_counter.size());
        sharedData.forbidden_move_counter[index * context.partition.k + m.to].fetch_add(1, std::memory_order_acq_rel);
      }
    }
    touched_edges.clear();
    move_edges_begin.clear();
  }

  template<typename FMStrategy>
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  bool LocalizedKWayFM<FMStrategy>::moveForbidden(PartitionedHypergraph& phg, Move& move) {
    if (sharedData.num_large_he == 0 || move.to == kInvalidPartition) {
      return false;
    }
    for (auto e : phg.incidentEdges(move.node)) {
      ASSERT(e < sharedData.num_edges_up_to.size());
      if (sharedData.num_edges_up_to[e + 1] - sharedData.num_edges_up_to[e] == 1) {
        ASSERT(phg.edgeSize(e) >= context.refinement.fm.large_he_threshold);
        size_t edge_id = sharedData.num_edges_up_to[e];
        ASSERT(edge_id < sharedData.num_large_he);
        ASSERT(edge_id * context.partition.k + move.to < sharedData.forbidden_move_counter.size());
        if (sharedData.forbidden_move_counter[edge_id * context.partition.k + move.to].load(std::memory_order_relaxed) >= context.refinement.fm.forbidden_move_theshold) {
          return true;
        }
      }
    }
    return false;
  }

  template<typename FMStrategy>
    template<typename PHG>
    MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
    void LocalizedKWayFM<FMStrategy>::syncMessageQueues(PHG &phg) {
      SearchID this_index = thisSearch - sharedData.nodeTracker.deactivatedNodeMarker - 1;
      size_t num_threads = context.shared_memory.num_threads;
      size_t mq_begin = this_index * num_threads;
      size_t mq_end = mq_begin + num_threads;
      Move m;
      for (size_t i = mq_begin; i < mq_end; ++i) {
        HypernodeID v;
        while (sharedData.messages[i].read(v)) {
          // use deduplicator to prevent uneeded pq updates
          if (neighborDeduplicator[v] != deduplicationTime && !sharedData.nodeTracker.isLocked(v)
              && sharedData.nodeTracker.searchOfNode[v] == thisSearch) {
            // this forces gain recalculation as we do not want to put moves into the mq
            m.from = sharedData.targetPart[v];
            fm_strategy.updateGain(phg, v, m);
            neighborDeduplicator[v] = deduplicationTime;
          }
        }
        fm_strategy.updatePQs(phg);
        /*sharedData.messages[i].clear();*/
      }
      updateNeighborDeduplicator();
      _local_moves_since_sync = 0;
    }

  template<typename FMStrategy>
  void LocalizedKWayFM<FMStrategy>::memoryConsumption(utils::MemoryTreeNode *parent) const {
    ASSERT(parent);

    utils::MemoryTreeNode *localized_fm_node = parent->addChild("Localized k-Way FM");

    utils::MemoryTreeNode *deduplicator_node = localized_fm_node->addChild("Deduplicator");
    deduplicator_node->updateSize(neighborDeduplicator.capacity() * sizeof(HypernodeID));
    utils::MemoryTreeNode *edges_to_activate_node = localized_fm_node->addChild("edgesWithGainChanges");
    edges_to_activate_node->updateSize(edgesWithGainChanges.capacity() * sizeof(HyperedgeID));

    utils::MemoryTreeNode *local_moves_node = parent->addChild("Local FM Moves");
    local_moves_node->updateSize(localMoves.capacity() * sizeof(std::pair<Move, MoveID>));

    fm_strategy.memoryConsumption(localized_fm_node);
    // TODO fm_strategy.memoryConsumptiom(..)
    /*
    utils::MemoryTreeNode* block_pq_node = localized_fm_node->addChild("Block PQ");
    block_pq_node->updateSize(blockPQ.size_in_bytes());
    utils::MemoryTreeNode* vertex_pq_node = localized_fm_node->addChild("Vertex PQ");
    for ( const VertexPriorityQueue& pq : vertexPQs ) {
      vertex_pq_node->updateSize(pq.size_in_bytes());
    }
     */

    deltaPhg.memoryConsumption(localized_fm_node);
  }

}   // namespace mt_kahypar


// instantiate templates
#include "mt-kahypar/partition/refinement/fm/strategies/gain_cache_strategy.h"
#include "mt-kahypar/partition/refinement/fm/strategies/gain_delta_strategy.h"
#include "mt-kahypar/partition/refinement/fm/strategies/recompute_gain_strategy.h"
#include <mt-kahypar/partition/refinement/fm/strategies/gain_cache_on_demand_strategy.h>

namespace mt_kahypar {
  template class LocalizedKWayFM<GainCacheStrategy>;
  template class LocalizedKWayFM<GainDeltaStrategy>;
  template class LocalizedKWayFM<RecomputeGainStrategy>;
  template class LocalizedKWayFM<GainCacheOnDemandStrategy>;
}
