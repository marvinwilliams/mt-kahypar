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
#include "mt-kahypar/partition/refinement/fm/scheduler_localized_kway_fm_core.h"
#include <queue>
namespace mt_kahypar {

  template<typename FMStrategy>
  class LocalSearchScheduler {
  public:
    explicit LocalSearchScheduler(const Context& context, const HypernodeID numNodes, FMSharedData& sharedData) :
      numNodes(numNodes),
      context(context),
      sharedData(sharedData),
      ets_fm([&] { return constructLocalizedKWayFMSearch(); }) { }

    void performLocalSearches(PartitionedHypergraph& phg, size_t numSeeds, size_t numSearches);

  private:

    SchedulerLocalizedKWayFM<FMStrategy> constructLocalizedKWayFMSearch() {
      return SchedulerLocalizedKWayFM<FMStrategy>(context,numNodes, sharedData);
    }

    void initSearches(PartitionedHypergraph& phg, size_t numSeeds, size_t numSearches);

    bool scheduleSearch();

  private:
    const HypernodeID numNodes;
    const Context& context;
    FMSharedData& sharedData;
    tbb::enumerable_thread_specific<SchedulerLocalizedKWayFM<FMStrategy>> ets_fm;

    std::priority_queue<std::pair<Gain, size_t>> local_searches;

    vec<SearchData<FMStrategy>> search_data;

    tbb::task_group tg;

    std::mutex m;
  };

}
