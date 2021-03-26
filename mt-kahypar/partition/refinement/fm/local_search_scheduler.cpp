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

#include "mt-kahypar/partition/refinement/fm/local_search_scheduler.h"

namespace mt_kahypar {
  template<typename FMStrategy>
  void LocalSearchScheduler<FMStrategy>::performLocalSearches(
    PartitionedHypergraph& phg, size_t numSeeds, size_t numSearches) {
      search_data = vec<SearchData<FMStrategy>>(numSearches);
      initSearches(phg, numSeeds, numSearches);
      while (!local_searches.empty()) {
        local_searches.pop();
      }
      tbb::task_group tg;
      for (size_t i = 0; i < search_data.size(); ++i) {
        Gain gain = search_data[i].gain;
        if (gain != invalidGain) {
          local_searches.emplace(gain, i);
        }
      }
      auto task = [&]() {
        auto& fm = ets_fm.local();
        while (sharedData.finishedTasks.load(std::memory_order_relaxed) < sharedData.finishedTasksLimit) {
          m.lock();
          if (local_searches.empty()) {
            m.unlock();
            break;
          }
          size_t search = local_searches.top().second;
          auto& data = search_data[search];
          local_searches.pop();
          m.unlock();
          data.scheduleStats.searches_started_and_resumed++;
          auto result = fm.resumeLocalSearch(phg, data);
          if (result.has_value()) { // reinsert to resume later
            m.lock();
            local_searches.emplace(result.value(), search);
            m.unlock();
          } else { // reinsert boundary vertices and reinsert search into pq
            sharedData.finishedTasks.fetch_add(1, std::memory_order_relaxed);
            fm.setup(phg, numSeeds, data);
            Gain gain = data.gain;
            if (gain != invalidGain) {
              m.lock();
              local_searches.emplace(gain, search);
              m.unlock();
            }
          }
        }
      };
      for (size_t i = 0; i < std::min(numSearches, context.shared_memory.num_threads); ++i) {
        tg.run(task);
      }
      tg.wait();
      if (debug) {
        ScheduleStats scheduleStats;
        for (auto& search : search_data) {
          search.scheduleStats.merge(scheduleStats);
        }
        LOG << scheduleStats.serialize();
        FMStats stats;
        collectStats(stats);
        LOG << stats.serialize();
      }

  }

  template<typename FMStrategy>
  void LocalSearchScheduler<FMStrategy>::collectStats(FMStats& stats) {
    for (auto& fm : ets_fm) {
      fm.stats.merge(stats);
    }
  }

  template<typename FMStrategy>
  void LocalSearchScheduler<FMStrategy>::initSearches(PartitionedHypergraph& phg, size_t numSeeds, size_t numSearches) {
    auto task = [&, numSeeds](const size_t search) {
      auto& fm = ets_fm.local();
      ASSERT(search_data[search].thisSearch == 0);
      fm.setup(phg, numSeeds, search_data[search]);
    };
    tbb::task_group tg;
    for (size_t i = 0; i < numSearches; ++i) {
      tg.run(std::bind(task, i));
    }
    tg.wait();
    /*tbb::parallel_for(size_t(0), numSearches, task);*/
  }
}

// instantiate templates
#include "mt-kahypar/partition/refinement/fm/strategies/gain_cache_strategy.h"
#include "mt-kahypar/partition/refinement/fm/strategies/gain_delta_strategy.h"
#include "mt-kahypar/partition/refinement/fm/strategies/recompute_gain_strategy.h"
#include <mt-kahypar/partition/refinement/fm/strategies/gain_cache_on_demand_strategy.h>

namespace mt_kahypar {
  template class LocalSearchScheduler<GainCacheStrategy>;
  template class LocalSearchScheduler<GainDeltaStrategy>;
  template class LocalSearchScheduler<RecomputeGainStrategy>;
  template class LocalSearchScheduler<GainCacheOnDemandStrategy>;
}
