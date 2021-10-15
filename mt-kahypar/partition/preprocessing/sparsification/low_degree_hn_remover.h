/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2016 Sebastian Schlag <sebastian.schlag@kit.edu>
 * Copyright (C) 2019 Tobias Heuer <tobias.heuer@kit.edu>
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

#include "tbb/enumerable_thread_specific.h"

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/parallel/stl/scalable_vector.h"

namespace mt_kahypar {

class LowDegreeHypernodeRemover {

 #define HIGH_DEGREE_HN_THRESHOLD ID(100)

 public:
  LowDegreeHypernodeRemover(const Context& context) :
    _context(context) { }

  LowDegreeHypernodeRemover(const LowDegreeHypernodeRemover&) = delete;
  LowDegreeHypernodeRemover & operator= (const LowDegreeHypernodeRemover &) = delete;

  LowDegreeHypernodeRemover(LowDegreeHypernodeRemover&&) = delete;
  LowDegreeHypernodeRemover & operator= (LowDegreeHypernodeRemover &&) = delete;

  HypernodeID removeLowDegreeHypernodes(Hypergraph& hg) {
    tbb::enumerable_thread_specific<HypernodeID> num_removed_low_degree_vertices(0);

    HypernodeID high_degree_threshold = highDegreeNodeThreshold(hg);
    hg.doParallelForAllNodes([&](const HypernodeID& hn) {
      if ( hg.nodeDegree(hn) < high_degree_threshold ) {
        bool removable = true;
        for ( const HyperedgeID& he : hg.incidentEdges(hn) ) {
          for ( const HypernodeID& pin : hg.pins(he) ) {
            if ( pin != hn && hg.nodeDegree(pin) < high_degree_threshold ) {
              removable = false;
              break;
            }
          }
          if ( !removable ) {
            break;
          }
        }

        if ( removable ) {
          ++num_removed_low_degree_vertices.local();
        }
      }
    });

    return num_removed_low_degree_vertices.combine(std::plus<HypernodeID>());
  }

  void restoreLowDegreeHypernodes(PartitionedHypergraph&) {

  }

  HypernodeID highDegreeNodeThreshold(const Hypergraph& hg) const {
    return std::max(static_cast<HypernodeID>(
      _context.preprocessing.high_degree_hn_threshold *
      hg.initialNumNodes()), HIGH_DEGREE_HN_THRESHOLD);
  }

 private:
  const Context& _context;
};

}  // namespace mt_kahypar
