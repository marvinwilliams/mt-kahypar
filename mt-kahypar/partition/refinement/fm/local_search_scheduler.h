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
#include <mt-kahypar/datastructures/priority_queue.h>
namespace mt_kahypar {

  template<typename FMStrategy>
  struct SearchData {
    SearchData(Context& context, HypernodeID numNodes, FMSharedData& sharedData) :
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
    explicit LocalSearchScheduler(const Context& context, HypernodeID numNodes, FMSharedData& sharedData) :
      ets_fm([&] { return constructLocalizedKWayFMSearch(context, numNodes, sharedData); }) {
        search_data.assign(context.shared_memory.num_threads + context.refinement.fm.num_additional_searches,
                           SearchData<FMStrategy>(context, numNodes, sharedData));
      }

    void performLocalSearches(PartitionedHypergraph& phg, size_t numSeeds);

    void initSearches();

  private:

    LocalizedKWayFM<FMStrategy> constructLocalizedKWayFMSearch(
      Context& context, HypernodeID numNodes, FMSharedData& sharedData) {
      return LocalizedKWayFM<FMStrategy>(context, numNodes, sharedData);
    }

    bool scheduleSearch();

  private:
    mt_kahypar::ds::Heap<Gain, size_t> local_searches;

    tbb::enumerable_thread_specific<LocalizedKWayFM<FMStrategy>> ets_fm;

    vec<SearchData<FMStrategy>> search_data;
  };

}
