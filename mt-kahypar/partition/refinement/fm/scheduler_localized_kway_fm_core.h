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


#pragma once

#include <mt-kahypar/partition/context.h>
#include <mt-kahypar/partition/metrics.h>

#include "mt-kahypar/datastructures/delta_partitioned_hypergraph.h"
#include "mt-kahypar/datastructures/sparse_map.h"
#include "mt-kahypar/partition/refinement/fm/fm_commons.h"
#include "mt-kahypar/partition/refinement/fm/stop_rule.h"

namespace mt_kahypar {


template<typename FMStrategy>
class SchedulerLocalizedKWayFM {
public:
  explicit SchedulerLocalizedKWayFM(const Context& context, HypernodeID numNodes, FMSharedData& sharedData) :
          context(context),
          k(context.partition.k),
          deltaPhg(context.partition.k),
          neighborDeduplicator(numNodes, 0),
          sharedData(sharedData)
          /*fm_strategy(context, numNodes, sharedData, {})*/
          { }

  bool setup(PartitionedHypergraph& phg, size_t numSeeds, SearchData<FMStrategy>& searchData);

  void memoryConsumption(utils::MemoryTreeNode* parent) const ;

  void checkDeltaMemory() {
    if (deltaPhg.combinedMemoryConsumption() > sharedData.deltaMemoryLimitPerThread) {
      sharedData.deltaExceededMemoryConstraints = true;
    }
  }

  /* TODO: make bool for return status aborted or can be resumed <07-03-21, @noahares> */
  std::optional<Gain> resumeLocalSearch(PartitionedHypergraph& phg, SearchData<FMStrategy>& search_data) {
    bool global_moves = context.refinement.fm.perform_moves_global || sharedData.deltaExceededMemoryConstraints;
    if (global_moves) {
      return internalFindMoves<false>(phg, search_data);
    } else {
      auto result = internalFindMoves<true>(phg, search_data);
      checkDeltaMemory();
      return result;
    }
  }

  FMStats stats;

private:

  // ! Performs localized FM local search on the delta partitioned hypergraph.
  // ! Moves made by this search are not immediately visible to other concurrent local searches.
  // ! The best prefix of moves is applied to the global partitioned hypergraph after the search finishes.
  //void internalFindMovesOnDeltaHypergraph(PartitionedHypergraph& phg, FMSharedData& sharedData);

  template<bool use_delta>
  std::optional<Gain> internalFindMoves(PartitionedHypergraph& phg, SearchData<FMStrategy>& search_data);


  template<typename PHG>
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  void acquireOrUpdateNeighbors(PHG& phg, const Move& move);


  // ! Makes moves applied on delta hypergraph visible on the global partitioned hypergraph.
  std::pair<Gain, size_t> applyBestLocalPrefixToSharedPartition(PartitionedHypergraph& phg,
                                                                const size_t best_index_locally_observed,
                                                                const Gain best_improvement_locally_observed,
                                                                bool apply_all_moves);

  // ! Rollback to the best improvement found during local search in case we applied moves
  // ! directly on the global partitioned hypergraph.
  void revertToBestLocalPrefix(PartitionedHypergraph& phg, size_t bestGainIndex);

  void applyMovesOntoDeltaPhg();

 private:

  const Context& context;

  // ! Number of blocks
  PartitionID k;

  // ! Wrapper around the global partitioned hypergraph, that allows
  // ! to perform moves non-visible for other local searches
  ds::DeltaPartitionedHypergraph<PartitionedHypergraph> deltaPhg;

  // ! Used after a move. Stores whether a neighbor of the just moved vertex has already been updated.
  vec<HypernodeID> neighborDeduplicator;
  HypernodeID deduplicationTime = 0;

  // ! Stores hyperedges whose pins's gains may have changed after vertex move
  vec<HyperedgeID> edgesWithGainChanges;

  FMSharedData& sharedData;

  /*FMStrategy fm_strategy;*/

  SearchData<FMStrategy>* searchData;

};

}
