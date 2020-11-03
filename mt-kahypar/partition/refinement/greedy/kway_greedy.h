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

#include "mt-kahypar/partition/context.h"

#include "mt-kahypar/datastructures/delta_partitioned_hypergraph.h"
#include "mt-kahypar/partition/refinement/fm/fm_commons.h"
#include "mt-kahypar/partition/refinement/fm/strategies/gain_cache_strategy.h"
#include "mt-kahypar/partition/refinement/i_refiner.h"
#include "mt-kahypar/partition/refinement/policies/gain_policy.h"

namespace mt_kahypar {

class KWayGreedy {

public:
  explicit KWayGreedy(const Context &c, HypernodeID numNodes,
                      FMSharedData &sharedData)
      : context(c), thisSearch(0), k(context.partition.k),
        neighborDeduplicator(numNodes, 0),
        fm_strategy(context, numNodes, sharedData, runStats),
        sharedData(sharedData), _gain(context) {}

  bool findMoves(PartitionedHypergraph &phg, size_t taskID);

  Gain getGain() { return _gain.delta(); }

  void memoryConsumption(utils::MemoryTreeNode *parent) const;

  FMStats stats;

private:
  void internalFindMoves(PartitionedHypergraph &phg);

  template <typename PHG>
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE void updateNeighbors(PHG &phg,
                                                          const Move &move);

  const Context &context;

  // ! Unique search id associated with the current local search
  SearchID thisSearch;

  // ! Number of blocks
  PartitionID k;

  // ! Local data members required for one localized search run
  // FMLocalData localData;
  vec<std::pair<Move, MoveID>> localMoves;

  // ! Used after a move. Stores whether a neighbor of the just moved vertex has
  // already been updated.
  vec<HypernodeID> neighborDeduplicator;
  HypernodeID deduplicationTime = 0;

  // ! Stores hyperedges whose pins's gains may have changed after vertex move
  vec<HyperedgeID> edgesWithGainChanges;

  FMStats runStats;

  GainCacheStrategy fm_strategy;

  FMSharedData &sharedData;

  GainPolicy<PartitionedHypergraph> _gain;
};

} // namespace mt_kahypar
