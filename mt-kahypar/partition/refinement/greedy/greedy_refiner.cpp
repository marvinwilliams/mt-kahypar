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

#include "mt-kahypar/partition/refinement/greedy/greedy_refiner.h"

namespace mt_kahypar {

  void BasicGreedyRefiner::initializeImpl(PartitionedHypergraph& phg) {
    // implement some initialization
  }

  /**
   *
   * @param phg The partition and hypergraph to refine
   * @param refinement_nodes Used for n-level uncoarsening. Ignore for now. It will be empty
   * @param metrics stores current km1 and imbalance
   * @return whether the refinement improved the objective function
   */
  bool BasicGreedyRefiner::refineImpl(PartitionedHypergraph& phg,
                                      const parallel::scalable_vector<HypernodeID>& refinement_nodes,
                                      kahypar::Metrics& metrics, double ) {

    // implement the refinement

    // don't forget to set the new imbalance and km1 values in the metrics object. you can ignore the cut value

    LOG << "You called greedy refinement";

    return false;
  }


}
