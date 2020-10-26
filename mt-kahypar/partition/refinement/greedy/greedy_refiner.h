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

#include "mt-kahypar/partition/refinement/i_refiner.h"


namespace mt_kahypar {

class BasicGreedyRefiner final : public IRefiner {

  static constexpr bool debug = false;
  static constexpr bool enable_heavy_assert = false;


public:

  BasicGreedyRefiner(const Hypergraph& hypergraph, const Context& c, const TaskGroupID taskGroupID)
  {
  }

  bool refineImpl(PartitionedHypergraph& phg,
                  const parallel::scalable_vector<HypernodeID>& refinement_nodes,
                  kahypar::Metrics& metrics,
                  double ) final ;

  void initializeImpl(PartitionedHypergraph& phg) final ;

};

} // namespace mt_kahypar
