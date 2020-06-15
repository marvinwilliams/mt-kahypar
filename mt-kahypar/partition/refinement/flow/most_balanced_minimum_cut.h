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
using KaHyParGraph = kahypar::ds::Graph;
using kahypar::ds::Edge;

template <typename FlowTypeTraits>
class MostBalancedMinimumCut {
 private:
  using FlowNetwork = typename FlowTypeTraits::FlowNetwork;
  using Scheduler = typename FlowTypeTraits::Scheduler;
  using Derived = typename FlowTypeTraits::MostBalancedMinimumCut;


 public:
  MostBalancedMinimumCut(const size_t initial_size) :
    _visited(initial_size),
    _graph_to_flow_network(initial_size, kinvalidFlowNetworkNode),
    _flow_network_to_graph(initial_size, kinvalidFlowNetworkNode),
    _scc_node_weight(initial_size, 0),
    _node_assignment(initial_size),
    _fixed_part_weights(std::make_pair(0,0)),
    _not_aquired_part_weights(std::make_pair(0,0)),
    _Q(),
    _sccs(initial_size) { }

  MostBalancedMinimumCut(const MostBalancedMinimumCut&) = delete;
  MostBalancedMinimumCut(MostBalancedMinimumCut&&) = delete;
  MostBalancedMinimumCut& operator= (const MostBalancedMinimumCut&) = delete;
  MostBalancedMinimumCut& operator= (MostBalancedMinimumCut&&) = delete;

  void mostBalancedMinimumCut(PartitionedHypergraph& hypergraph, FlowNetwork& flow_network,
                              const Context& context,
                              const PartitionID block_0, const PartitionID block_1, Scheduler& scheduler, double& new_imbalance) {
    reset();

    // Mark all reachable nodes from source and sink set as invalid
    markAllReachableNodesAsVisited<true>(hypergraph, flow_network);
    markAllReachableNodesAsVisited<false>(hypergraph, flow_network);

    //check if fixed weights already exceed the limit
    HypernodeWeight not_aq_first = scheduler.get_not_aquired_weight(block_0, block_1);
    HypernodeWeight not_aq_second = scheduler.get_not_aquired_weight(block_1, block_0);
    std::pair<HypernodeWeight, HypernodeWeight> fixed_per_part = std::make_pair(_fixed_part_weights.first + not_aq_first,
    _fixed_part_weights.second + not_aq_second);
    if(fixed_per_part.first > context.partition.max_part_weights[block_0] ||
     fixed_per_part.second > context.partition.max_part_weights[block_1]){
       new_imbalance = 1;
       return;
     }

    // Build residual graph
    KaHyParGraph residual_graph(std::move(buildResidualGraph(hypergraph, flow_network)));

    // Find strongly connected components
    findStronglyConnectedComponents(residual_graph);

    // Contract strongly connected components
    auto contraction = residual_graph.contractClusters();
    KaHyParGraph dag(std::move(contraction.first));
    std::vector<NodeID> contraction_mapping = std::move(contraction.second);

    // Build mapping from contracted graph to flow network
    parallel::scalable_vector<parallel::scalable_vector<NodeID> > scc_to_flow_network(dag.numNodes(), parallel::scalable_vector<NodeID>());
    for (const NodeID& u_og : residual_graph.nodes()) {
      const NodeID flow_u_og = _graph_to_flow_network.get(u_og);
      if (flow_network.isHypernode(flow_u_og)) {
        scc_to_flow_network[contraction_mapping[u_og]].push_back(flow_u_og);
        _scc_node_weight.update(contraction_mapping[u_og], hypergraph.nodeWeight(flow_u_og));
      }
    }

    // Calculate in degrees of nodes in DAG graph for topological ordering
    parallel::scalable_vector<size_t> in_degree(dag.numNodes(), 0);
    for (const NodeID& u : dag.nodes()) {
      for (const Edge& e : dag.incidentEdges(u)) {
        const NodeID v = e.target_node;
        if (u != v) {
          in_degree[v]++;
        }
      }
    }

    // Find most balanced minimum cut
    parallel::scalable_vector<NodeID> topological_order(dag.numNodes(), 0);
    parallel::scalable_vector<bool> best_partition_id(dag.numNodes(), false);

    std::pair<HypernodeWeight, HypernodeWeight> aquired_part_weight = get_aquired_part_weight(block_0, block_1, scheduler);
    HypernodeWeight initial_current_first = aquired_part_weight.first + aquired_part_weight.second - _fixed_part_weights.second;
    std::pair<HypernodeWeight, HypernodeWeight> initial_current_part_weight = std::make_pair(initial_current_first, _fixed_part_weights.second);

    not_aq_first = scheduler.get_not_aquired_weight(block_0, block_1);
    not_aq_second = scheduler.get_not_aquired_weight(block_1, block_0);
    _not_aquired_part_weights = std::make_pair(not_aq_first, not_aq_second);

    double best_imbalance = imbalance(context, initial_current_part_weight);

    for (size_t i = 0; i < 20; ++i) {
      // Compute random topological order
      topologicalSort(dag, in_degree, topological_order);

      // Sweep through topological order and find best imbalance
      parallel::scalable_vector<bool> tmp_partition_id(dag.numNodes(), false);
      double tmp_best_imbalance = best_imbalance;//metrics::localBlockImbalance(hypergraph, context, block_0, block_1);

      std::pair<HypernodeWeight, HypernodeWeight> current_part_weight = initial_current_part_weight;

      for (size_t idx = 0; idx < topological_order.size(); ++idx) {
        const NodeID u = topological_order[idx];
        tmp_partition_id[u] = true;
        current_part_weight.first -= _scc_node_weight.get(u);
        current_part_weight.second += _scc_node_weight.get(u);
        double cur_imbalance = static_cast<Derived*>(this)->imbalance(context, current_part_weight);

        if (cur_imbalance > tmp_best_imbalance) {
          tmp_partition_id[u] = false;
          break;
        }
        tmp_best_imbalance = cur_imbalance;
      }

      if (tmp_best_imbalance < best_imbalance) {
        best_imbalance = tmp_best_imbalance;
        best_partition_id = tmp_partition_id;
        ASSERT(current_part_weight.first + current_part_weight.second == aquired_part_weight.first + aquired_part_weight.second);
      }
    }

    DBG << "Best imbalance: " << best_imbalance;

    new_imbalance = best_imbalance;

    // add most balanced minimum cut to assignment
    for (const NodeID& u : dag.nodes()) {
      for (const NodeID& v : scc_to_flow_network[u]) {
        _node_assignment.set(v, best_partition_id[u]);
      }
    }

     ASSERT([&]() {
      kahypar::ds::FastResetFlagArray<> moved(hypergraph.initialNumNodes());
      for(const HypernodeID& hn : flow_network.hypernodes()) {
          const PartitionID from = hypergraph.partID(hn);
          const PartitionID to = _node_assignment[hn] ? block_1: block_0;
          if (from != to) {
              hypergraph.changeNodePart(hn, from, to);
              moved.set(hn, true);
          }
      }

      for (const HypernodeID& hn : flow_network.hypernodes()) {
          if(moved[hn]){
            const PartitionID from = hypergraph.partID(hn);
            const PartitionID to = from == block_0? block_1: block_0;
            hypergraph.changeNodePart(hn, from, to);
            moved.set(hn, true);
          }
      }
      return true;
    } (), "Most balanced minimum cut failed!");

  }

  kahypar::ds::FastResetFlagArray<>& get_assignment(){
    return _node_assignment;
  }


  protected:
  std::pair<HypernodeWeight, HypernodeWeight> get_aquired_part_weight(size_t block_0,size_t block_1, Scheduler& scheduler){
    return scheduler.get_aquired_part_weight(block_0, block_1);
  }

 private:
  static constexpr bool debug = false;

  void reset() {
    _visited.reset();
    _graph_to_flow_network.resetUsedEntries();
    _flow_network_to_graph.resetUsedEntries();
    _scc_node_weight.resetUsedEntries();
    _node_assignment.reset();
    _fixed_part_weights = std::make_pair(0, 0);
    _not_aquired_part_weights = std::make_pair(0, 0);
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
  void markAllReachableNodesAsVisited(PartitionedHypergraph& hypergraph, FlowNetwork& flow_network) {
    auto start_set_iterator = sourceSet ? flow_network.sources() : flow_network.sinks();
    for (const NodeID& node : start_set_iterator) {
      _Q.push(node);
      _visited.set(node, true);
      assign_and_add_weight(hypergraph, node, sourceSet, flow_network);
    }

    while (!_Q.empty()) {
      const NodeID u = _Q.front();

      _Q.pop();

      if (flow_network.interpreteHyperedge(u, sourceSet)) {
        const HyperedgeID he = flow_network.mapToHyperedgeID(u);
        for (const HypernodeID& pin : hypergraph.pins(he)) {
          if (flow_network.containsHypernode(pin)) {
            if (flow_network.isRemovedHypernode(pin)) {
              if(!_visited[pin]){
                _visited.set(pin, true);
                if (!sourceSet){
                  _node_assignment.set(pin, true);
                  _fixed_part_weights.second += hypergraph.nodeWeight(pin);
                }else{
                  _fixed_part_weights.first += hypergraph.nodeWeight(pin);
                }
              }
            }
          }
        }
      }

      for(ds::FlowEdge & e:flow_network.incidentEdges(u)){
        const ds::FlowEdge& reverse_edge = flow_network.reverseEdge(e);
        const NodeID v = e.target;
        if (!_visited[v]) {
          if ((sourceSet && flow_network.residualCapacity(e)) ||
              (!sourceSet && flow_network.residualCapacity(reverse_edge))) {
            _Q.push(v);
            _visited.set(v, true);
            assign_and_add_weight(hypergraph, v, sourceSet, flow_network);
          }
        }
      }
    }
  }


  void assign_and_add_weight(PartitionedHypergraph& hypergraph, HypernodeID node, bool sourceSet, FlowNetwork& flow_network){
    if (flow_network.interpreteHypernode(node)){
      if (!sourceSet){
        _node_assignment.set(node, true);
        _fixed_part_weights.second += hypergraph.nodeWeight(node);
      } else {
        _fixed_part_weights.first += hypergraph.nodeWeight(node);
      }
    }
  }

  KaHyParGraph buildResidualGraph(PartitionedHypergraph& hypergraph, FlowNetwork& flow_network) {
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

    parallel::scalable_vector<parallel::scalable_vector<Edge> > adj_list(cur_graph_node, parallel::scalable_vector<Edge>());
    for (const NodeID& node : flow_network.nodes()) {
      if (!_visited[node]) {
        const NodeID source = _flow_network_to_graph.get(node);

        for(ds::FlowEdge & flow_edge:flow_network.incidentEdges(node)){
          const NodeID target = flow_edge.target;
          if (flow_network.residualCapacity(flow_edge) && !_visited[target]) {
            Edge e;
            e.target_node = _flow_network_to_graph.get(target);
            e.weight = 1.0;
            adj_list[source].push_back(e);
          }
        }
      }
    }

    for (const HypernodeID& hn : flow_network.removedHypernodes()) {
      if (!_visited[hn]) {
        const NodeID hn_node = _flow_network_to_graph.get(hn);
        for (const HyperedgeID& he : hypergraph.incidentEdges(hn)) {
          const NodeID in_he = _flow_network_to_graph.get(flow_network.mapToIncommingHyperedgeID(he));
          const NodeID out_he = _flow_network_to_graph.get(flow_network.mapToOutgoingHyperedgeID(he));
          if (in_he != kinvalidFlowNetworkNode) {
            Edge e;
            e.target_node = in_he;
            e.weight = 1.0;
            adj_list[hn_node].push_back(e);
          }
          if (out_he != kinvalidFlowNetworkNode) {
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

    return KaHyParGraph(adj_array, edges);
  }

  void findStronglyConnectedComponents(KaHyParGraph& g) {
    _sccs.compute(g);
  }


  void topologicalSort(const KaHyParGraph& g,
                       parallel::scalable_vector<size_t> in_degree,
                       parallel::scalable_vector<NodeID>& topological_order) {
    parallel::scalable_vector<NodeID> start_nodes;
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

  inline double imbalance(const Context& context, const std::pair<HypernodeWeight, HypernodeWeight>& part_weight) {
    const HypernodeWeight weight_part0 = part_weight.first;
    const HypernodeWeight weight_part1 = part_weight.second;
    const double imbalance_part0 = ((_not_aquired_part_weights.first + weight_part0) /
                                    static_cast<double>(context.partition.perfect_balance_part_weights[0]));
    const double imbalance_part1 = ((_not_aquired_part_weights.second + weight_part1) /
                                    static_cast<double>(context.partition.perfect_balance_part_weights[1]));
    return std::max(imbalance_part0, imbalance_part1) - 1.0;
  }

  kahypar::ds::FastResetFlagArray<> _visited;
  kahypar::ds::FastResetArray<NodeID> _graph_to_flow_network;
  kahypar::ds::FastResetArray<NodeID> _flow_network_to_graph;
  kahypar::ds::FastResetArray<HypernodeWeight> _scc_node_weight;
  //false indicates block_0, true block_1
  kahypar::ds::FastResetFlagArray<> _node_assignment;
  std::pair<HypernodeWeight, HypernodeWeight> _fixed_part_weights;
protected:
  std::pair<HypernodeWeight, HypernodeWeight> _not_aquired_part_weights;
private:
  std::queue<NodeID> _Q;  // BFS queue
  StronglyConnectedComponents _sccs;
};

template <typename FlowTypeTraits>
class MatchingMostBalancedMinimumCut : public MostBalancedMinimumCut<FlowTypeTraits> {
  using Scheduler = typename FlowTypeTraits::Scheduler;
  using FlowNetwork = typename FlowTypeTraits::FlowNetwork;
  using Base = MostBalancedMinimumCut<FlowTypeTraits>;

  using Base::Base;
};

template <typename FlowTypeTraits>
class OptMostBalancedMinimumCut : public MostBalancedMinimumCut<FlowTypeTraits> {
  using Scheduler = typename FlowTypeTraits::Scheduler;
  using FlowNetwork = typename FlowTypeTraits::FlowNetwork;
  using Base = MostBalancedMinimumCut<FlowTypeTraits>;

  using Base::Base;
};

}  // namespace mt_kahypar
