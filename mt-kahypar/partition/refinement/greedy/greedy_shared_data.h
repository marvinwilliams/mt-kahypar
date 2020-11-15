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

#include "mt-kahypar/datastructures/hypergraph_common.h"
#include "mt-kahypar/parallel/hold_barrier.h"
#include <condition_variable>
#include <mutex>
namespace mt_kahypar {

using HypernodeIDMessageMatrix = vec<vec<HypernodeID>>;
struct GreedySharedData {
  HypernodeIDMessageMatrix messages;
  parallel::HoldBarrier hold_barrier;
  GreedySharedData(size_t num_threads)
      : messages(num_threads * num_threads, vec<HypernodeID>()),
        hold_barrier(num_threads) {}
};
} // namespace mt_kahypar
