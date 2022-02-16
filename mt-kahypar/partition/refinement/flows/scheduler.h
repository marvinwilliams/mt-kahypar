/*******************************************************************************
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2021 Tobias Heuer <tobias.heuer@kit.edu>
 *
 * Mt-KaHyPar is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Mt-KaHyPar is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Mt-KaHyPar.  If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/

#pragma once

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/partition/refinement/i_refiner.h"
#include "mt-kahypar/partition/refinement/flows/quotient_graph.h"
#include "mt-kahypar/partition/refinement/flows/refiner_adapter.h"
#include "mt-kahypar/partition/refinement/flows/problem_construction.h"
#include "mt-kahypar/parallel/atomic_wrapper.h"
#include "mt-kahypar/utils/refinement_stats.h"

namespace mt_kahypar {

class FlowRefinementScheduler final : public IRefiner {

  static constexpr bool debug = false;
  static constexpr bool enable_heavy_assert = false;

  struct PartWeightUpdateResult {
    bool is_balanced = true;
    PartitionID overloaded_block = kInvalidPartition;
    HypernodeWeight overload_weight = 0;
  };

public:
  explicit FlowRefinementScheduler(const Hypergraph& hg,
                                   const Context& context) :
    _phg(nullptr),
    _context(context),
    _quotient_graph(hg, context),
    _refiner(hg, context),
    _constructor(hg, context),
    _was_moved(hg.initialNumNodes(), uint8_t(false)),
    _part_weights_lock(),
    _part_weights(context.partition.k, 0),
    _max_part_weights(context.partition.k, 0),
    _apply_moves_lock(),
    _flow_stats(nullptr) { }

  FlowRefinementScheduler(const FlowRefinementScheduler&) = delete;
  FlowRefinementScheduler(FlowRefinementScheduler&&) = delete;

  FlowRefinementScheduler & operator= (const FlowRefinementScheduler &) = delete;
  FlowRefinementScheduler & operator= (FlowRefinementScheduler &&) = delete;

  /**
   * Applies the sequence of vertex moves to the partitioned hypergraph.
   * The method ensures that the move sequence does not violate
   * the balance constaint and not worsen solution quality.
   * Returns, improvement in solution quality.
   */
  HyperedgeWeight applyMoves(const SearchID search_id,
                             MoveSequence& sequence);

  /**
   * Returns the current weight of each block.
   * Note, we do not want that the underlying refiner (ILP and Flows)
   * see partially updated part weight information. Thus, we perform
   * part weight updates for a move sequence as a transaction, which
   * we protect with a spin lock.
   */
  vec<HypernodeWeight> partWeights() {
    _part_weights_lock.lock();
    vec<HypernodeWeight> _copy_part_weights(_part_weights);
    _part_weights_lock.unlock();
    return _copy_part_weights;
  }

private:
  bool refineImpl(PartitionedHypergraph& phg,
                  const vec<HypernodeID>& refinement_nodes,
                  Metrics& metrics,
                  double time_limit) final;

  void initializeImpl(PartitionedHypergraph& phg) final;

  PartWeightUpdateResult partWeightUpdate(const vec<HypernodeWeight>& part_weight_deltas,
                                          const bool rollback);

  std::string blocksOfSearch(const SearchID search_id) {
    const BlockPair blocks = _quotient_graph.getBlockPair(search_id);
    return "(" + std::to_string(blocks.i) + "," + std::to_string(blocks.j) + ")";
  }

  PartitionedHypergraph* _phg;
  const Context& _context;

  // ! Contains information of all cut hyperedges between the
  // ! blocks of the partition
  QuotientGraph _quotient_graph;

  // ! Maintains the flow refiner instances
  FlowRefinerAdapter _refiner;

  // ! Responsible for construction of an flow problems
  ProblemConstruction _constructor;

  // ! For each vertex it store wheather the corresponding vertex
  // ! was moved or not
  vec<uint8_t> _was_moved;

  // ! Maintains the part weights of each block
  SpinLock _part_weights_lock;
  vec<HypernodeWeight> _part_weights;
  vec<HypernodeWeight> _max_part_weights;

  SpinLock _apply_moves_lock;

  vec<utils::FlowStats>* _flow_stats;
};

}  // namespace kahypar
