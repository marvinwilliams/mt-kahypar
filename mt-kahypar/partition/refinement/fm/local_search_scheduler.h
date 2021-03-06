/*******************************************************************************
 * This file is part of MT-KaHyPar.
 *
 * Copyright (C) 2021 Noah Wahl <noah.ares@yahoo.de>
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
#include "mt-kahypar/partition/refinement/fm/fm_commons.h"
#include "mt-kahypar/partition/refinement/fm/localized_kway_fm_core.h"
#include <queue>
namespace mt_kahypar {

  template<typename FMStrategy>
  struct SearchData {
    // each thread gets a pointer to the data of its current search
    SearchData(const Context& context, const HypernodeID numNodes, const FMSharedData& sharedData) :
      thisSearch(0),
      fm_strategy(context, numNodes, sharedData, runStats)
    { }
    SearchID thisSearch;
    vec<std::pair<Move, MoveID>> localMoves;
    FMStrategy fm_strategy; // TODO: share this or sperate per search? For now each search gets one
    FMStats runStats;
  };

  template<typename FMStrategy>
  class LocalSearchScheduler {
  public:
    explicit LocalSearchScheduler(const Context& context, const HypernodeID numNodes, const FMSharedData& sharedData) :
      numNodes(numNodes),
      context(context),
      sharedData(sharedData),
      ets_fm([&] { return constructLocalizedKWayFMSearch(); }) { }

    void performLocalSearches(PartitionedHypergraph& phg, size_t numSeeds, size_t numSearches);

  private:

    LocalizedKWayFM<FMStrategy> constructLocalizedKWayFMSearch() {
      return LocalizedKWayFM<FMStrategy>(context,numNodes, sharedData);
    }

    void initSearches(PartitionedHypergraph& phg, size_t numSeeds);

    bool scheduleSearch();

  private:
    const HypernodeID numNodes;
    const Context& context;
    const FMSharedData& sharedData;
    tbb::enumerable_thread_specific<LocalizedKWayFM<FMStrategy>> ets_fm;

    std::priority_queue<std::pair<Gain, size_t>> local_searches;

    vec<SearchData<FMStrategy>> search_data;
  };

}
