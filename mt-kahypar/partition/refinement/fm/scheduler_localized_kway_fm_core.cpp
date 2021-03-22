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

#include "mt-kahypar/partition/refinement/fm/scheduler_localized_kway_fm_core.h"
#include "mt-kahypar/partition/refinement/fm/strategies/km1_gains.h"

namespace mt_kahypar {

  template<typename FMStrategy>
  bool SchedulerLocalizedKWayFM<FMStrategy>::setup(PartitionedHypergraph& phg, size_t numSeeds, SearchData<FMStrategy>& _searchData) {
    searchData = &_searchData;
    fm_strategy.setRunStats(searchData->runStats);
    searchData->localMoves.clear();
    searchData->nodes.clear();
    if (searchData->thisSearch == 0) {
      searchData->thisSearch = ++sharedData.nodeTracker.highestActiveSearchID;
    }
    ASSERT(searchData->thisSearch - sharedData.nodeTracker.deactivatedNodeMarker <= context.shared_memory.num_threads + context.refinement.fm.additional_searches_factor);

    auto seeds = sharedData.shared_refinement_nodes.try_pop(numSeeds);
    Gain max_gain = invalidGain;
    if (seeds) {
      Km1GainComputer gc(sharedData.numParts);
      for (HypernodeID u : *seeds) {
        if (sharedData.nodeTracker.tryAcquireNode(u, searchData->thisSearch)) {
          searchData->nodes.push_back(u);
          // REVIEW can look in gain cache? ask strategy for the gain value
          auto [target, gain] = gc.computeBestTargetBlock(phg, u, context.partition.max_part_weights);
          max_gain = std::max(max_gain, gain);
        }
      }
    }

    searchData->gain = max_gain;
    if (searchData->nodes.size() > 0) {
      if (sharedData.deltaExceededMemoryConstraints) {
        deltaPhg.dropMemory();
      }
      return true;
    } else {
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
  void SchedulerLocalizedKWayFM<FMStrategy>::acquireOrUpdateNeighbors(PHG& phg, const Move& move) {
    // Note: In theory we should acquire/update all neighbors. It just turned out that this works fine
    // Actually: only vertices incident to edges with gain changes can become new boundary vertices.
    // Vertices that already were boundary vertices, can still be considered later since they are in the task queue
    // --> actually not that bad
    for (HyperedgeID e : edgesWithGainChanges) {
      if (phg.edgeSize(e) < context.partition.ignore_hyperedge_size_threshold) {
        for (HypernodeID v : phg.pins(e)) {
          if (neighborDeduplicator[v] != deduplicationTime) {
            SearchID searchOfV = sharedData.nodeTracker.searchOfNode[v].load(std::memory_order_acq_rel);
            if (searchOfV == searchData->thisSearch) {
              fm_strategy.updateGain(phg, v, move);
            } else if (sharedData.nodeTracker.tryAcquireNode(v, searchData->thisSearch)) {
              fm_strategy.insertIntoPQ(phg, v, searchOfV);
            }
            neighborDeduplicator[v] = deduplicationTime;
          }
        }
      }
    }
    edgesWithGainChanges.clear();

    if (++deduplicationTime == 0) {
      neighborDeduplicator.assign(neighborDeduplicator.size(), 0);
      deduplicationTime = 1;
    }
  }


  template<typename FMStrategy>
  template<bool use_delta>
  std::optional<Gain> SchedulerLocalizedKWayFM<FMStrategy>::internalFindMoves(PartitionedHypergraph& phg, SearchData<FMStrategy>& _searchData) {
    StopRule stopRule(phg.initialNumNodes());
    Move move;
    searchData = &_searchData;
    setFMStrategy(phg);
    if constexpr (use_delta) {
      deltaPhg.clear();
      deltaPhg.setPartitionedHypergraph(&phg);
      applyMovesOntoDeltaPhg();
    }

    auto delta_func = [&](const HyperedgeID he,
                          const HyperedgeWeight edge_weight,
                          const HypernodeID,
                          const HypernodeID pin_count_in_from_part_after,
                          const HypernodeID pin_count_in_to_part_after) {
      // Gains of the pins of a hyperedge can only change in the following situations.
      if (pin_count_in_from_part_after == 0 || pin_count_in_from_part_after == 1 ||
          pin_count_in_to_part_after == 1 || pin_count_in_to_part_after == 2) {
        edgesWithGainChanges.push_back(he);
      }

      if constexpr (use_delta) {
        fm_strategy.deltaGainUpdates(deltaPhg, he, edge_weight, move.from, pin_count_in_from_part_after,
                                     move.to, pin_count_in_to_part_after);
      } else {
        fm_strategy.deltaGainUpdates(phg, he, edge_weight, move.from, pin_count_in_from_part_after,
                                     move.to, pin_count_in_to_part_after);
      }

    };

    // we can almost make this function take a generic partitioned hypergraph
    // we would have to add the success func to the interface of DeltaPhg (and then ignore it there...)
    // and do the local rollback outside this function


    HypernodeWeight heaviestPartWeight = 0;
    HypernodeWeight fromWeight = 0, toWeight = 0;
    Gain next_gain = 0;

    while (!stopRule.searchShouldStop()
           && sharedData.finishedTasks.load(std::memory_order_relaxed) < sharedData.finishedTasksLimit) {

      Gain last_gain = next_gain;
      next_gain = fm_strategy.getNextMoveGain(phg);
      if (next_gain == kInvalidGain) break;
      // Idea: if more than half of the scheduled moves have been performed and the next move would be the first nagative gain move, reschedule
      bool preemptive_reschedule = last_gain >= 0 && next_gain < 0
                                   && searchData->num_moves > (context.refinement.fm.max_moves_before_reschedule / 2);
      if (preemptive_reschedule || searchData->num_moves > context.refinement.fm.max_moves_before_reschedule) {
        searchData->num_moves = 0;
        searchData->scheduleStats.reschedules++;
        if (preemptive_reschedule) {
          searchData->scheduleStats.preemptive_reschedules++;
        }
        fm_strategy.resetPQs(searchData->nodes);
        return next_gain;
      }

      if constexpr (use_delta) {
        if (!fm_strategy.findNextMove(deltaPhg, move)) break;
      } else {
        if (!fm_strategy.findNextMove(phg, move)) break;
      }
      ASSERT(!sharedData.nodeTracker.isLocked(move.node));

      sharedData.nodeTracker.deactivateNode(move.node, searchData->thisSearch);
      MoveID move_id = std::numeric_limits<MoveID>::max();
      bool moved = false;
      if (move.to != kInvalidPartition) {
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
        searchData->runStats.moves++;
        searchData->num_moves++;
        searchData->estimatedImprovement += move.gain;
        searchData->localMoves.emplace_back(move, move_id);
        stopRule.update(move.gain);
        const bool improved_km1 = searchData->estimatedImprovement > searchData->bestImprovement;
        const bool improved_balance_less_equal_km1 = searchData->estimatedImprovement >= searchData->bestImprovement
                                                     && fromWeight == heaviestPartWeight
                                                     && toWeight + phg.nodeWeight(move.node) < heaviestPartWeight;

        if (improved_km1 || improved_balance_less_equal_km1) {
          stopRule.reset();
          searchData->bestImprovement = searchData->estimatedImprovement;
          searchData->bestImprovementIndex = searchData->localMoves.size();

          if constexpr (use_delta) {
            applyBestLocalPrefixToSharedPartition(phg, searchData->bestImprovementIndex, searchData->bestImprovement, true /* apply all moves */);
            searchData->bestImprovementIndex = 0;
            searchData->localMoves.clear();
            deltaPhg.clear();   // clear hashtables, save memory :)
          }
        }

        if constexpr (use_delta) {
          acquireOrUpdateNeighbors(deltaPhg, move);
        } else {
          acquireOrUpdateNeighbors(phg, move);
        }
      }

    }

    if constexpr (use_delta) {
      std::tie(searchData->bestImprovement, searchData->bestImprovementIndex) =
              applyBestLocalPrefixToSharedPartition(phg, searchData->bestImprovementIndex, searchData->bestImprovement, false);
    } else {
      revertToBestLocalPrefix(phg, searchData->bestImprovementIndex);
    }

    searchData->runStats.estimated_improvement = searchData->bestImprovement;
    fm_strategy.clearPQs(searchData->bestImprovementIndex);
    searchData->reset();
    searchData->runStats.merge(stats);
    return {};
  }


  template<typename FMStrategy>
  std::pair<Gain, size_t> SchedulerLocalizedKWayFM<FMStrategy>::applyBestLocalPrefixToSharedPartition(
          PartitionedHypergraph& phg,
          const size_t best_index_locally_observed,
          const Gain best_improvement_locally_observed,
          bool apply_all_moves) {

    Gain improvement_from_attributed_gains = 0;
    Gain attributed_gain = 0;

    auto delta_gain_func = [&](const HyperedgeID he,
                               const HyperedgeWeight edge_weight,
                               const HypernodeID edge_size,
                               const HypernodeID pin_count_in_from_part_after,
                               const HypernodeID pin_count_in_to_part_after) {
      attributed_gain += km1Delta(he, edge_weight, edge_size,
                                  pin_count_in_from_part_after, pin_count_in_to_part_after);
    };

    // Apply move sequence to original hypergraph and update gain values
    Gain best_improvement_from_attributed_gains = 0;
    size_t best_index_from_attributed_gains = 0;
    for (size_t i = 0; i < best_index_locally_observed; ++i) {
      assert(i < searchData->localMoves.size());
      Move& local_move = searchData->localMoves[i].first;
      MoveID& move_id = searchData->localMoves[i].second;
      attributed_gain = 0;

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

      attributed_gain = -attributed_gain; // delta func yields negative sum of improvements, i.e. negative values mean improvements
      improvement_from_attributed_gains += attributed_gain;
      ASSERT(move_id != std::numeric_limits<MoveID>::max());
      if (improvement_from_attributed_gains >= best_improvement_from_attributed_gains) {
        best_improvement_from_attributed_gains = improvement_from_attributed_gains;
        best_index_from_attributed_gains = i;
      }
    }

    searchData->runStats.local_reverts += searchData->localMoves.size() - best_index_locally_observed;
    if (!apply_all_moves && best_index_from_attributed_gains != best_index_locally_observed) {
      searchData->runStats.best_prefix_mismatch++;
    }

    // kind of double rollback, if attributed gains say we overall made things worse
    if (!apply_all_moves && improvement_from_attributed_gains < 0) {
      // always using the if-branch gave similar results
      searchData->runStats.local_reverts += best_index_locally_observed - best_index_from_attributed_gains + 1;
      for (size_t i = best_index_from_attributed_gains + 1; i < best_index_locally_observed; ++i) {
        Move& m = sharedData.moveTracker.getMove(searchData->localMoves[i].second);

        if constexpr (FMStrategy::uses_gain_cache) {
          phg.changeNodePartWithGainCacheUpdate(m.node, m.to, m.from);
        } else {
          phg.changeNodePart(m.node, m.to, m.from);
        }

        m.invalidate();
      }
      return std::make_pair(best_improvement_from_attributed_gains, best_index_from_attributed_gains);
    } else {
      return std::make_pair(best_improvement_locally_observed, best_index_locally_observed);
    }
  }

  template<typename FMStrategy>
  void SchedulerLocalizedKWayFM<FMStrategy>::revertToBestLocalPrefix(PartitionedHypergraph& phg, size_t bestGainIndex) {
    searchData->runStats.local_reverts += searchData->localMoves.size() - bestGainIndex;
    while (searchData->localMoves.size() > bestGainIndex) {
      Move& m = sharedData.moveTracker.getMove(searchData->localMoves.back().second);
      if constexpr (FMStrategy::uses_gain_cache) {
        phg.changeNodePartWithGainCacheUpdate(m.node, m.to, m.from);
      } else {
        phg.changeNodePart(m.node, m.to, m.from);
      }
      m.invalidate();
      searchData->localMoves.pop_back();
    }
  }

  template<typename FMStrategy>
  void SchedulerLocalizedKWayFM<FMStrategy>::applyMovesOntoDeltaPhg() {
    for (auto& m : searchData->localMoves) {
      Move& move = m.first;
      /* TODO: is this right now? <22-03-21, @noahares> */
      // REVIEW nope. call the fm_strategy.deltaGainUpdates function in the delta_func lambda for changeNodePart
      auto delta_func = [&](const HyperedgeID he, const HyperedgeWeight edge_weight, const HypernodeID edge_size,
            const HypernodeID pin_count_in_from_part_after, const HypernodeID pin_count_in_to_part_after) {
        fm_strategy.deltaGainUpdates(
          deltaPhg, he, edge_weight, move.from, pin_count_in_from_part_after, move.to, pin_count_in_to_part_after);
      };

      deltaPhg.changeNodePart(
        move.node, move.from, move.to, context.partition.max_part_weights[move.to], delta_func);
    }
  }

  /* TODO: dont call insert multiple times, but use some sort of heapify <18-03-21, @noahares> */
  template<typename FMStrategy>
  void SchedulerLocalizedKWayFM<FMStrategy>::setFMStrategy(PartitionedHypergraph& phg) {
    fm_strategy.setRunStats(searchData->runStats);
    fm_strategy.insertAll(phg, searchData->nodes);
  }

  template<typename FMStrategy>
  void SchedulerLocalizedKWayFM<FMStrategy>::memoryConsumption(utils::MemoryTreeNode *parent) const {
    ASSERT(parent);

    utils::MemoryTreeNode *localized_fm_node = parent->addChild("Localized k-Way FM");

    utils::MemoryTreeNode *deduplicator_node = localized_fm_node->addChild("Deduplicator");
    deduplicator_node->updateSize(neighborDeduplicator.capacity() * sizeof(HypernodeID));
    utils::MemoryTreeNode *edges_to_activate_node = localized_fm_node->addChild("edgesWithGainChanges");
    edges_to_activate_node->updateSize(edgesWithGainChanges.capacity() * sizeof(HyperedgeID));

    utils::MemoryTreeNode *local_moves_node = parent->addChild("Local FM Moves");
    local_moves_node->updateSize(searchData->localMoves.capacity() * sizeof(std::pair<Move, MoveID>));

    fm_strategy.memoryConsumption(localized_fm_node);
    deltaPhg.memoryConsumption(localized_fm_node);
  }

}   // namespace mt_kahypar


// instantiate templates
#include "mt-kahypar/partition/refinement/fm/strategies/gain_cache_strategy.h"
#include "mt-kahypar/partition/refinement/fm/strategies/gain_delta_strategy.h"
#include "mt-kahypar/partition/refinement/fm/strategies/recompute_gain_strategy.h"
#include <mt-kahypar/partition/refinement/fm/strategies/gain_cache_on_demand_strategy.h>

namespace mt_kahypar {
  template class SchedulerLocalizedKWayFM<GainCacheStrategy>;
  template class SchedulerLocalizedKWayFM<GainDeltaStrategy>;
  template class SchedulerLocalizedKWayFM<RecomputeGainStrategy>;
  template class SchedulerLocalizedKWayFM<GainCacheOnDemandStrategy>;
}
