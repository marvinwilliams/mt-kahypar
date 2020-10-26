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

#include "kahypar/datastructure/kway_priority_queue.h"
#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/partition/refinement/i_refiner.h"
#include <mt-kahypar/datastructures/concurrent_bucket_map.h>
#include <mt-kahypar/datastructures/priority_queue.h>
#include <mt-kahypar/parallel/work_stack.h>

namespace mt_kahypar {

class BasicGreedyRefiner final : public IRefiner {

  static constexpr bool debug = false;
  static constexpr bool enable_heavy_assert = false;

  using KWayRefinementPQ =
      kahypar::ds::KWayPriorityQueue<HypernodeID, Gain,
                                     std::numeric_limits<Gain>>;

  /**
   * A hyperedge can be in three states during FM refinement: FREE, LOOSE and
   * LOCKED. Initial all hyperedges are FREE. Once we move a vertex incident to
   * the hyperedge, the hyperedge becomes LOOSE. If we move an other vertex in
   * the opposite direction the hyperedge becomes LOCKED. LOCKED hyperedges have
   * the property that they can not be removed from cut and we can therefore
   * skip delta gain updates.
   */
  enum HEState {
    FREE = std::numeric_limits<PartitionID>::max() - 1,
    LOCKED = std::numeric_limits<PartitionID>::max(),
  };

  /**
   * INACTIVE = Initial state of a vertex
   * ACTIVE = Vertex is a border node and inserted into the PQ
   * MOVED = Vertex was moved during local search
   */
  enum class VertexState { INACTIVE, ACTIVE, MOVED };

  /* TODO: constructor with init <21-10-20, @noahares> */
  struct BasicGreedySharedData {

    // ! Nodes to initialize the localized FM searches with
    WorkContainer<HypernodeID> refinement_nodes;

    // ! PQ handles shared by all threads (each vertex is only held by one
    // thread)
    vec<PosT> vertex_pq_handles;

    // ! num parts
    PartitionID num_parts;

    // ! Stores the designated target part of a vertex, i.e. the part with the
    // highest gain to which moving is feasible
    vec<PartitionID> target_part;

    // bucket map for storing hypernodes by their maximum gain
    ds::ConcurrentBucketMap<HypernodeID> gain_buckets;

    // stores the highest gain bucket that contains a vertex
    Gain best_gain;

    // ! Stop parallel refinement if finishedTasks > finishedTasksLimit to
    // avoid long-running single searches
    CAtomic<size_t> finished_tasks;
    size_t finished_tasks_limit = std::numeric_limits<size_t>::max();

    bool release_nodes = true;

    BasicGreedySharedData(size_t num_nodes, const Context &context)
        : refinement_nodes(), vertex_pq_handles(),
          num_parts(context.partition.k), target_part(), gain_buckets() {}
  };

public:
  BasicGreedyRefiner(const Hypergraph &hypergraph, const Context &c,
                     const TaskGroupID taskGroupID)
      : shared_data(hypergraph.initialNumNodes(), c), _hg(hypergraph),
        _initial_num_nodes(hypergraph.initialNumNodes()), _context(c),
        _taskGroupID(taskGroupID),
        _border_vertices(hypergraph.initialNumNodes()), _pq(c.partition.k),
        _vertex_state(hypergraph.initialNumNodes(), VertexState::INACTIVE),
        _he_state(hypergraph.initialNumEdges(), HEState::FREE) {}

  bool
  refineImpl(PartitionedHypergraph &phg,
             const parallel::scalable_vector<HypernodeID> &refinement_nodes,
             kahypar::Metrics &metrics, double) final;

  void initializeImpl(PartitionedHypergraph &phg) final;

  void initBorderVertices(
      PartitionedHypergraph &phg,
      const parallel::scalable_vector<HypernodeID> &refinement_nodess);

  BasicGreedySharedData shared_data;

private:
  bool _is_initialized = false;
  const Hypergraph &_hg;
  const HypernodeID _initial_num_nodes;
  const Context &_context;
  const TaskGroupID _taskGroupID;
  parallel::scalable_vector<HypernodeID> _nodes;
  BorderVertexTracker _border_vertices;
  /* TODO: do I need a pq? <23-10-20, @noahares> */
  KWayRefinementPQ _pq;
  parallel::scalable_vector<VertexState> _vertex_state;
  parallel::scalable_vector<PartitionID> _he_state;

  /*! \brief Find the best partition to move a vertex to
   *
   * \param phg the partitioned hypergraph
   * \param u the vertex to move
   * \return PartitionID, Gain pair of best partition to move u to with the gain
   */
  std::pair<PartitionID, Gain>
  findBestToPartition(const PartitionedHypergraph &phg, const HypernodeID u);

  /*! \brief Update Neighbor gains after moving a node
   *
   * \param phg the partitioned hypergraph
   * \param u the node that was moved
   */
  void updateNeighbors(PartitionedHypergraph &phg, HypernodeID u);
};

} // namespace mt_kahypar
