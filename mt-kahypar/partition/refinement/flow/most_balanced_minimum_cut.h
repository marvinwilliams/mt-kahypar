/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2018 Sebastian Schlag <sebastian.schlag@kit.edu>
 * Copyright (C) 2018 Tobias Heuer <tobias.heuer@live.com>
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

#include <algorithm>
#include <array>
#include <queue>
#include <stack>
#include <string>
#include <utility>
#include <vector>

#include "external_tools/kahypar/kahypar/datastructure/fast_reset_array.h"
#include "external_tools/kahypar/kahypar/datastructure/fast_reset_flag_array.h"
#include "external_tools/kahypar/kahypar/datastructure/sparse_set.h"
#include "external_tools/kahypar/kahypar/datastructure/graph.h"

#include "mt-kahypar/datastructures/flow_network.h"


#include "mt-kahypar/definitions.h"
#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/partition/metrics.h"
#include "mt-kahypar/partition/refinement/flow/strongly_connected_components.h"
#include "mt-kahypar/utils/randomize.h"

namespace mt_kahypar {
using kahypar::ds::Graph;
using kahypar::ds::Edge;

template <typename TypeTraits, typename FlowTypeTraits>
class MostBalancedMinimumCut {
 private:
  using HyperGraph = typename TypeTraits::template PartitionedHyperGraph<>;
  using FlowNetwork = typename FlowTypeTraits::FlowNetwork;
  using Scheduler = typename FlowTypeTraits::Scheduler;
  using Derived = typename FlowTypeTraits::MostBalancedMinimumCut;


 public:
  MostBalancedMinimumCut(const size_t initial_size) :
    _visited(initial_size),
    _graph_to_flow_network(initial_size, FlowNetwork::kInvalidNode),
    _flow_network_to_graph(initial_size, FlowNetwork::kInvalidNode),
    _scc_node_weight(initial_size, 0),
    _Q(),
    _sccs(initial_size) { }

  MostBalancedMinimumCut(const MostBalancedMinimumCut&) = delete;
  MostBalancedMinimumCut(MostBalancedMinimumCut&&) = delete;
  MostBalancedMinimumCut& operator= (const MostBalancedMinimumCut&) = delete;
  MostBalancedMinimumCut& operator= (MostBalancedMinimumCut&&) = delete;

  void mostBalancedMinimumCut(HyperGraph& hypergraph, FlowNetwork& flow_network,
                              const Context& context,
                              const PartitionID block_0, const PartitionID block_1, Scheduler& scheduler) {
    reset();

    // Mark all reachable nodes from source and sink set as invalid
    markAllReachableNodesAsVisited<true>(hypergraph, flow_network, block_0, block_1);
    markAllReachableNodesAsVisited<false>(hypergraph, flow_network, block_0, block_1);

    // Build residual graph
    Graph residual_graph(std::move(buildResidualGraph(hypergraph, flow_network)));

    // Find strongly connected components
    findStronglyConnectedComponents(residual_graph);

    // Contract strongly connected components
    auto contraction = residual_graph.contractClusters();
    Graph dag(std::move(contraction.first));
    std::vector<NodeID> contraction_mapping = std::move(contraction.second);

    // Build mapping from contracted graph to flow network
    std::vector<std::vector<NodeID> > scc_to_flow_network(dag.numNodes(), std::vector<NodeID>());
    for (const NodeID& u_og : residual_graph.nodes()) {
      const NodeID flow_u_og = _graph_to_flow_network.get(u_og);
      if (flow_network.isHypernode(flow_u_og)) {
        scc_to_flow_network[contraction_mapping[u_og]].push_back(flow_u_og);
        _scc_node_weight.update(contraction_mapping[u_og], hypergraph.nodeWeight(hypergraph.globalNodeID(flow_u_og)));
      }
    }

    // Calculate in degrees of nodes in DAG graph for topological ordering
    std::vector<size_t> in_degree(dag.numNodes(), 0);
    for (const NodeID& u : dag.nodes()) {
      for (const Edge& e : dag.incidentEdges(u)) {
        const NodeID v = e.target_node;
        if (u != v) {
          in_degree[v]++;
        }
      }
    }

    // Find most balanced minimum cut
    std::vector<NodeID> topological_order(dag.numNodes(), 0);
    std::vector<PartitionID> best_partition_id(dag.numNodes(), block_0);
    std::vector<HypernodeWeight> aquired_part_weight = flow_network.get_aquired_part_weight(hypergraph, block_0, block_1);
    double best_imbalance = Derived::imbalance(hypergraph, context, aquired_part_weight, scheduler, block_0, block_1);

    DBG << "Start Most Balanced Minimum Cut (Bipartition = {" << block_0 << "," << block_1 << "}";
    DBG << "Initial imbalance: " << V(Derived::imbalance(hypergraph, context, aquired_part_weight, scheduler, block_0, block_1));

    for (size_t i = 0; i < 20; ++i) {
      // Compute random topological order
      topologicalSort(dag, in_degree, topological_order);

      // Sweep through topological order and find best imbalance
      std::vector<PartitionID> tmp_partition_id(dag.numNodes(), block_0);
      double tmp_best_imbalance = best_imbalance;//metrics::localBlockImbalance(hypergraph, context, block_0, block_1);

      std::vector<HypernodeWeight> part_weight = aquired_part_weight;

      for (size_t idx = 0; idx < topological_order.size(); ++idx) {
        const NodeID u = topological_order[idx];
        tmp_partition_id[u] = block_1;
        part_weight[0] -= _scc_node_weight.get(u);
        part_weight[1] += _scc_node_weight.get(u);
        double cur_imbalance = Derived::imbalance(hypergraph, context, part_weight, scheduler, block_0, block_1);

        if (cur_imbalance > tmp_best_imbalance) {
          tmp_partition_id[u] = block_0;
          break;
        }
        tmp_best_imbalance = cur_imbalance;
      }

      if (tmp_best_imbalance < best_imbalance) {
        best_imbalance = tmp_best_imbalance;
        best_partition_id = tmp_partition_id;
      }
    }

    DBG << "Best imbalance: " << best_imbalance;

    ASSERT([&]() {
        std::vector<HypernodeWeight> assert_aquired_part_weight = flow_network.get_aquired_part_weight(hypergraph, block_0, block_1);
        const double imbalance_before = Derived::imbalance(hypergraph, context, aquired_part_weight, scheduler, block_0, block_1);

        std::vector<NodeID> topological_order(dag.numNodes(), 0);
        std::vector<NodeID> part_before(dag.numNodes(), block_0);
        topologicalSort(dag, in_degree, topological_order);
        for (const NodeID& u_og : topological_order) {
          for (const NodeID& v_og : scc_to_flow_network[u_og]) {
            const PartitionID from = hypergraph.partID(hypergraph.globalNodeID(v_og));
            const PartitionID to = best_partition_id[u_og];
            if (from != to) {
              bool success = hypergraph.changeNodePart(hypergraph.globalNodeID(v_og), from, to);
              unused(success);
              ASSERT(success);
              part_before[u_og] = from;
            }
          }
        }

        // Rollback hypernode assignment
        for (const NodeID& u : dag.nodes()) {
          for (const NodeID& v_og : scc_to_flow_network[u]) {
            const PartitionID from = hypergraph.partID(hypergraph.globalNodeID(v_og));
            const PartitionID to = part_before[u];
            if (from != to) {
              bool success = hypergraph.changeNodePart(hypergraph.globalNodeID(v_og), from, to);
              unused(success);
              ASSERT(success);
            }
          }
        }

        //const HyperedgeWeight metric = metrics::objective(hypergraph, context.partition.objective);
        //metric != metric_before ||
        if (Derived::imbalance(hypergraph, context, aquired_part_weight, scheduler, block_0, block_1) != imbalance_before) {
          LOG << "Restoring original partition failed!";
          LOG << V(imbalance_before) << V(Derived::imbalance(hypergraph, context, aquired_part_weight, scheduler, block_0, block_1));
          return false;
        }

        return true;
      } (), "Most balanced minimum cut failed!");

    // Assign most balanced minimum cut
    for (const NodeID& u : dag.nodes()) {
      for (const NodeID& v_og : scc_to_flow_network[u]) {
        const PartitionID from = hypergraph.partID(hypergraph.globalNodeID(v_og));
        const PartitionID to = best_partition_id[u];
        if (from != to) {
          bool success = hypergraph.changeNodePart(hypergraph.globalNodeID(v_og), from, to);
          unused(success);
          ASSERT(success);
        }
      }
    }

    /*ASSERT(best_imbalance == metrics::localBlockImbalance(hypergraph, context, block_0, block_1),
           "Best imbalance didn't match with current imbalance"
           << V(best_imbalance) << V(metrics::localBlockImbalance(hypergraph, context, block_0, block_1)));*/
  }


  protected:
  std::vector<HypernodeWeight> get_aquired_part_weight(HyperGraph& hypergraph, size_t block_0,size_t block_1, FlowNetwork& flow_network){
    return static_cast<Derived*>(this)->get_aquired_part_weight_impl(hypergraph,block_0, block_1, flow_network);
  }

 private:
  static constexpr bool debug = false;

  void reset() {
    _visited.reset();
    _graph_to_flow_network.resetUsedEntries();
    _flow_network_to_graph.resetUsedEntries();
    _scc_node_weight.resetUsedEntries();
  }

  /**
   * Executes a BFS starting from the source (sourceSet = true)
   * or sink set (sourceSet = false). Touched nodes by BFS
   * are marked as visited and are not considered during
   * residual graph building step.
   *
   * @t_param sourceSet Indicates, if BFS start from source or sink set
   */
  template <bool sourceSet>
  void markAllReachableNodesAsVisited(HyperGraph& hypergraph, FlowNetwork& flow_network,
                                      const PartitionID block_0, const PartitionID block_1) {
    auto start_set_iterator = sourceSet ? flow_network.sources() : flow_network.sinks();
    for (const NodeID& node : start_set_iterator) {
      _Q.push(node);
      _visited.set(node, true);
    }

    while (!_Q.empty()) {
      const NodeID u_og = _Q.front();

      _Q.pop();

      if (flow_network.interpreteHypernode(u_og)) {
        const HypernodeID u = hypergraph.globalNodeID(u_og);
        if (!sourceSet) {
          const PartitionID from = hypergraph.partID(u);
          if (from == block_0) {
            bool success = hypergraph.changeNodePart(u, block_0, block_1);
            unused(success);
            ASSERT(success);
          }
        }
      } else if (flow_network.interpreteHyperedge(u_og, sourceSet)) {
        const HyperedgeID he_og = flow_network.mapToHyperedgeID(u_og);
        const HyperedgeID he = hypergraph.globalEdgeID(he_og);
        for (const HypernodeID& pin : hypergraph.pins(he)) {
          if (flow_network.containsHypernode(hypergraph, pin)) {
            if (!sourceSet) {
              PartitionID from = hypergraph.partID(pin);
              if (from == block_0) {
                bool success = hypergraph.changeNodePart(pin, block_0, block_1);
                unused(success);
                ASSERT(success);
              }
            }
            if (flow_network.isRemovedHypernode(hypergraph, pin)) {
              _visited.set(hypergraph.originalNodeID(pin), true);
            }
          }
        }
      }

      size_t next_idx = flow_network.getFirstFlowEdge(u_og);
      
      while(next_idx != FlowNetwork::kInvalidNode){
        FlowEdge & e = flow_network.getEdge(next_idx);
        const FlowEdge& reverse_edge = flow_network.reverseEdge(e);
        const NodeID v_og = e.target;
        if (!_visited[v_og]) {
          if ((sourceSet && flow_network.residualCapacity(e)) ||
              (!sourceSet && flow_network.residualCapacity(reverse_edge))) {
            _Q.push(v_og);
            _visited.set(v_og, true);
          }
        }
        next_idx = e.nextEdge;
      }
    }
  }

  Graph buildResidualGraph(HyperGraph& hypergraph, FlowNetwork& flow_network) {
    size_t cur_graph_node = 0;
    for (const NodeID& node : flow_network.nodes()) {
      if (!_visited[node]) {
        _graph_to_flow_network.set(cur_graph_node, node);
        _flow_network_to_graph.set(node, cur_graph_node++);
      }
    }

    for (const HypernodeID& hn : flow_network.removedHypernodes()) {
      if (!_visited[hn]) {
        _graph_to_flow_network.set(cur_graph_node, hn);
        _flow_network_to_graph.set(hn, cur_graph_node++);
      }
    }

    std::vector<std::vector<Edge> > adj_list(cur_graph_node, std::vector<Edge>());
    for (const NodeID& node : flow_network.nodes()) {
      if (!_visited[node]) {
        const NodeID source = _flow_network_to_graph.get(node);

        size_t next_idx = flow_network.getFirstFlowEdge(node);
      
        while(next_idx != FlowNetwork::kInvalidNode){
          FlowEdge & flow_edge = flow_network.getEdge(next_idx);
          const NodeID target = flow_edge.target;
          if (flow_network.residualCapacity(flow_edge) && !_visited[target]) {
            Edge e;
            e.target_node = _flow_network_to_graph.get(target);
            e.weight = 1.0;
            adj_list[source].push_back(e);
          }
          next_idx = flow_edge.nextEdge;
        }
      }
    }

    for (const HypernodeID& hn_og : flow_network.removedHypernodes()) {
      if (!_visited[hn_og]) {
        const NodeID hn_node = _flow_network_to_graph.get(hn_og);
        for (const HyperedgeID& he : hypergraph.incidentEdges(hypergraph.globalNodeID(hn_og))) {
          const NodeID in_he = _flow_network_to_graph.get(flow_network.mapToIncommingHyperedgeID(hypergraph, he));
          const NodeID out_he = _flow_network_to_graph.get(flow_network.mapToOutgoingHyperedgeID(hypergraph, he));
          if (in_he != FlowNetwork::kInvalidNode) {
            Edge e;
            e.target_node = in_he;
            e.weight = 1.0;
            adj_list[hn_node].push_back(e);
          }
          if (out_he != FlowNetwork::kInvalidNode) {
            Edge e;
            e.target_node = hn_node;
            e.weight = 1.0;
            adj_list[out_he].push_back(e);
          }
        }
      }
    }

    std::vector<NodeID> adj_array(1, 0);
    std::vector<Edge> edges;
    for (NodeID u = 0; u < cur_graph_node; ++u) {
      for (const Edge& e : adj_list[u]) {
        edges.push_back(e);
      }
      adj_array.push_back(adj_array[u] + adj_list[u].size());
    }

    return Graph(adj_array, edges);
  }

  void findStronglyConnectedComponents(Graph& g) {
    _sccs.compute(g);
  }


  void topologicalSort(const Graph& g,
                       std::vector<size_t> in_degree,
                       std::vector<NodeID>& topological_order) {
    std::vector<NodeID> start_nodes;
    for (const NodeID& u : g.nodes()) {
      if (in_degree[u] == 0) {
        start_nodes.push_back(u);
      }
    }
    utils::Randomize::instance().shuffleVector(start_nodes);
    for (const NodeID& u : start_nodes) {
      _Q.push(u);
    }

    size_t idx = 0;
    while (!_Q.empty()) {
      const NodeID u = _Q.front();
      _Q.pop();
      topological_order[idx++] = u;
      for (const Edge& e : g.incidentEdges(u)) {
        const NodeID v = e.target_node;
        if (u != v) {
          in_degree[v]--;
          if (in_degree[v] == 0) {
            _Q.push(v);
          }
        }
      }
    }

    ASSERT(idx == g.numNodes(), "Topological sort failed!" << V(idx) << V(g.numNodes()));
  }

  kahypar::ds::FastResetFlagArray<> _visited;
  kahypar::ds::FastResetArray<NodeID> _graph_to_flow_network;
  kahypar::ds::FastResetArray<NodeID> _flow_network_to_graph;
  kahypar::ds::FastResetArray<HypernodeWeight> _scc_node_weight;

  std::queue<NodeID> _Q;  // BFS queue
  StronglyConnectedComponents _sccs;
};

template <typename TypeTraits, typename FlowTypeTraits>
class MatchingMostBalancedMinimumCut : public MostBalancedMinimumCut<TypeTraits, FlowTypeTraits> {
  using HyperGraph = typename TypeTraits::template PartitionedHyperGraph<>;
  using Scheduler = typename FlowTypeTraits::Scheduler;
  using FlowNetwork = typename FlowTypeTraits::FlowNetwork;
  using Base = MostBalancedMinimumCut<TypeTraits, FlowTypeTraits>;

  using Base::Base;

  public:
  std::vector<HypernodeWeight> get_aquired_part_weight_impl(HyperGraph& hypergraph, size_t block_0,size_t block_1, FlowNetwork& flow_network){
    return std::vector<HypernodeWeight> (2, 0);
  }

  static inline double imbalance(HyperGraph& hypergraph, const Context& context, const std::vector<HypernodeWeight>& part_weight,
   Scheduler & scheduler, size_t block_0, size_t block_1) {
    unused(scheduler);
    unused(part_weight);
    return metrics::localBlockImbalance(hypergraph, context, block_0, block_1);
  }
};

template <typename TypeTraits, typename FlowTypeTraits>
class OptMostBalancedMinimumCut : public MostBalancedMinimumCut<TypeTraits, FlowTypeTraits> {
  using HyperGraph = typename TypeTraits::template PartitionedHyperGraph<>;
  using Scheduler = typename FlowTypeTraits::Scheduler;
  using FlowNetwork = typename FlowTypeTraits::FlowNetwork;
  using Base = MostBalancedMinimumCut<TypeTraits, FlowTypeTraits>;

  using Base::Base;

  public:
  std::vector<HypernodeWeight> get_aquired_part_weight_impl(HyperGraph& hypergraph, size_t block_0,size_t block_1, FlowNetwork& flow_network){
    return flow_network.get_aquired_part_weight(hypergraph,block_0, block_1);
  }

  static inline double imbalance(HyperGraph& hypergraph, const Context& context, const std::vector<HypernodeWeight>& part_weight,
  Scheduler & scheduler, size_t block_0, size_t block_1) {
    unused(hypergraph);
    const HypernodeWeight weight_part0 = part_weight[0];
    const HypernodeWeight weight_part1 = part_weight[1];
    const double imbalance_part0 = ((scheduler.get_not_aquired_weight(block_0, block_1) + weight_part0) /
                                    static_cast<double>(context.partition.perfect_balance_part_weights[0]));
    const double imbalance_part1 = ((scheduler.get_not_aquired_weight(block_1, block_0) + weight_part1) /
                                    static_cast<double>(context.partition.perfect_balance_part_weights[1]));
    return std::max(imbalance_part0, imbalance_part1) - 1.0;
  }
};

}  // namespace mt_kahypar
