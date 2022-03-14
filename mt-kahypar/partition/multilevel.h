/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2020 Tobias Heuer <tobias.heuer@kit.edu>
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

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/partition/context.h"

namespace mt_kahypar::multilevel {

// ! Performs multilevel partitioning on the given hypergraph
// ! in TBB blocking-style.
PartitionedHypergraph partition(Hypergraph& hypergraph,
                                              const Context& context,
                                              const bool top_level,
                                              const bool vcycle = false);
// ! Performs multilevel partitioning on the given hypergraph
// ! in TBB continuation-style.
// ! Note, the final partitioned hypergraph is moved into the
// ! passed partitioned hypergraph object.
void partition_async(Hypergraph& hypergraph, PartitionedHypergraph& partitioned_hypergraph,
                                   const Context& context,
                                   const bool top_level,
                                   tbb::task* parent);


}  // namespace mt_kahypar
