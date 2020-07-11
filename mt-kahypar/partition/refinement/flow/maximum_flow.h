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
#include <queue>
#include <unordered_map>
#include <utility>
#include <vector>

#include "external_maximum_flow/ibfs/ibfs.h"

#include "kahypar/datastructure/fast_reset_array.h"
#include "kahypar/datastructure/fast_reset_flag_array.h"
#include "kahypar/datastructure/sparse_set.h"
#include "kahypar/meta/mandatory.h"
#include "kahypar/meta/typelist.h"

#include "mt-kahypar/datastructures/flow_network.h"
#include "mt-kahypar/definitions.h"
#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/partition/metrics.h"
#include "mt-kahypar/utils/randomize.h"
#include "mt-kahypar/partition/refinement/flow/most_balanced_minimum_cut.h"

namespace mt_kahypar {
template <typename FlowTypeTraits>
class MaximumFlow {
using FlowNetwork = typename FlowTypeTraits::FlowNetwork;
using Scheduler = typename FlowTypeTraits::Scheduler;
using MostBalancedMinimumCut = typename FlowTypeTraits::MostBalancedMinimumCut;

 public:
  MaximumFlow(const HypernodeID initial_size, const HypernodeID initial_num_nodes) :
    _mbmc(initial_size),
    _new_imbalance(0),
    _original_part_id(initial_num_nodes, 0) { }

  virtual ~MaximumFlow() { }

  MaximumFlow(const MaximumFlow&) = delete;
  MaximumFlow& operator= (const MaximumFlow&) = delete;

  MaximumFlow(MaximumFlow&&) = delete;
  MaximumFlow& operator= (MaximumFlow&&) = delete;


  virtual Flow maximumFlow(PartitionedHypergraph& hypergraph, FlowNetwork& flow_network) = 0;

  HyperedgeWeight minimumSTCut(PartitionedHypergraph& hypergraph, FlowNetwork& flow_network,
                               const Context& context,
                               const PartitionID block_0, const PartitionID block_1, Scheduler & scheduler,
                               const HyperedgeWeight cut_flow_network_before) {
    if (flow_network.isTrivialFlow()) {
      return kInfty;
    }

    utils::Timer::instance().start_timer("maxFlow", "Calculating max Flow ", true);
    const HyperedgeWeight cut = maximumFlow(hypergraph, flow_network);
    utils::Timer::instance().stop_timer("maxFlow");

    if(cut < cut_flow_network_before){
      utils::Timer::instance().start_timer("MBMC", "Finding MBMC ", true);
      _mbmc.mostBalancedMinimumCut(hypergraph, flow_network, context, block_0, block_1, scheduler, _new_imbalance);
      utils::Timer::instance().stop_timer("MBMC");
    }
    return cut;
  }

  double get_new_imbalance(){
    return _new_imbalance;
  }

  kahypar::ds::FastResetFlagArray<>& get_assignment(){
    return _mbmc.get_assignment();
  }

  protected:
  MostBalancedMinimumCut _mbmc;
  double _new_imbalance;

  parallel::scalable_vector<PartitionID> _original_part_id;
};


template <typename FlowTypeTraits>
class IBFS : public MaximumFlow<FlowTypeTraits>{
  using Base = MaximumFlow<FlowTypeTraits>;
  using FlowGraph = maxflow::IBFSGraph;
 private:
  using FlowNetwork = typename FlowTypeTraits::FlowNetwork;


 public:
  IBFS(const HypernodeID initial_size, const HypernodeID initial_num_nodes) :
    Base(initial_size, initial_num_nodes),
    _flow_graph(FlowGraph::IB_INIT_FAST),
    _flow_network_mapping(initial_size, 0) { }

  ~IBFS() = default;

  IBFS(const IBFS&) = delete;
  IBFS& operator= (const IBFS&) = delete;

  IBFS(IBFS&&) = delete;
  IBFS& operator= (IBFS&&) = delete;

  Flow maximumFlow(PartitionedHypergraph& hypergraph, FlowNetwork& flow_network) {
    unused(hypergraph);
    mapToExternalFlowNetwork(flow_network);

    _flow_graph.computeMaxFlow();
    const Flow max_flow = _flow_graph.getFlow();

    //assign flow to flow network
    FlowGraph::Arc* a = _flow_graph.arcs;
    while (a != _flow_graph.arcEnd) {
      const Flow flow = a->flowEdge->capacity - a->rCap;
      if (flow != 0) {
        a->flowEdge->increaseFlow(flow);
      }
      a++;
    }
    return max_flow;
  }

 private:
  template <typename T>
  FRIEND_TEST(AMaximumFlow, ChecksIfAugmentingPathExist);
  template <typename T>
  FRIEND_TEST(AMaximumFlow, AugmentAlongPath);

  void mapToExternalFlowNetwork(FlowNetwork& flow_network) {
    _flow_graph.initSize(flow_network.numNodes(), flow_network.numEdges() - flow_network.numUndirectedEdges());
    const Flow infty = flow_network.totalWeightHyperedges();
    NodeID cur_id = 0;

    for (const NodeID& node : flow_network.nodes()) {
      _flow_graph.addNode(cur_id,
                          flow_network.isIdSource(node) ? infty : 0,
                          flow_network.isIdSink(node) ? infty : 0);
      _flow_network_mapping[node] = cur_id;
      cur_id++;
    }

    //add Edges
    bool reverse = false;
    for (ds::FlowEdge& edge : flow_network.edges()) {
      if(reverse){
        reverse = false;
      }else{
        const NodeID u = _flow_network_mapping[edge.source];
        const NodeID v = _flow_network_mapping[edge.target];
        const Capacity c = edge.capacity;
        ds::FlowEdge& rev_edge = flow_network.reverseEdge(edge);
        const Capacity rev_c = rev_edge.capacity;
        _flow_graph.addEdge(u, v, c, rev_c, &edge, &rev_edge);
        reverse = true;
      }
    }

    _flow_graph.initGraph();
  }

  FlowGraph _flow_graph;
  parallel::scalable_vector<NodeID> _flow_network_mapping;
};
}  // namespace mt_kahypar
