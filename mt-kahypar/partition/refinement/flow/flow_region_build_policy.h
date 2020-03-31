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
#include <vector>

#include "kahypar/datastructure/fast_reset_flag_array.h"
#include "mt-kahypar/datastructures/flow_network.h"
#include "kahypar/meta/policy_registry.h"
#include "kahypar/meta/typelist.h"
#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/utils/randomize.h"

namespace mt_kahypar {
template< typename TypeTraits>
class FlowRegionBuildPolicy : public kahypar::meta::PolicyBase {

  private:
    using HyperGraph = typename TypeTraits::template PartitionedHyperGraph<>;


 public:
  template <class Network = Mandatory, typename Scheduler>
  static inline HypernodeID bfs(HyperGraph& hg,
                                Network& flow_network,
                                std::vector<HypernodeID>& start_nodes,
                                const PartitionID part,
                                const PartitionID other_part,
                                const HypernodeWeight max_part_weight,
                                kahypar::ds::FastResetFlagArray<>& visited,
                                Scheduler& scheduler) {
    visited.reset();
    utils::Randomize::instance().shuffleVector(start_nodes);
    std::queue<HypernodeID> Q;
    HypernodeWeight queue_weight = 0;
    for (const HypernodeID& hn : start_nodes) {
      if (queue_weight + hg.nodeWeight(hn) <= max_part_weight) {
        Q.push(hn);
        queue_weight += hg.nodeWeight(hn);
        visited.set(hg.originalNodeID(hn), true);
      }else{
        scheduler.releaseNode(hn);
      }
    }

    HypernodeID num_hypernodes_added = 0;
    const size_t num_hypernodes = hg.initialNumNodes();
    while (!Q.empty()) {
      const HypernodeID hn = Q.front();
      Q.pop();

      flow_network.addHypernode(hg, hn);
      ++num_hypernodes_added;

      for (const HyperedgeID& he : hg.incidentEdges(hn)) {
        if (!visited[num_hypernodes + hg.originalEdgeID(he)]) {
          for (const HypernodeID& pin : hg.pins(he)) {
            if (!visited[hg.originalNodeID(pin)] && hg.partID(pin) == part &&
                queue_weight + hg.nodeWeight(pin) <= max_part_weight) {
              if(scheduler.tryAquireNode(pin)){
                Q.push(pin);
                queue_weight += hg.nodeWeight(pin);
              }
              visited.set(hg.originalNodeID(pin), true);
            }
          }
          visited.set(num_hypernodes + hg.originalEdgeID(he), true);
        }
      }
    }
    scheduler.aquire_block_weight(part, other_part, queue_weight);
    return num_hypernodes_added;
  }
};

template< typename TypeTraits>
class CutBuildPolicy : public FlowRegionBuildPolicy<TypeTraits> {
  private:
    using HyperGraph = typename TypeTraits::template PartitionedHyperGraph<>;
 public:
  template <class Network = Mandatory, typename Scheduler>
  inline static void buildFlowNetwork(HyperGraph& hg,
                                      const Context& context,
                                      Network& flow_network,
                                      std::vector<HyperedgeID>& cut_hes,
                                      const double alpha,
                                      const PartitionID block_0,
                                      const PartitionID block_1,
                                      kahypar::ds::FastResetFlagArray<>& visited,
                                      Scheduler & scheduler) {
    visited.reset();
    std::vector<HypernodeID> start_nodes_block_0;
    std::vector<HypernodeID> start_nodes_block_1;
    for (const HyperedgeID he : cut_hes) {
      //TODO: why does this fail sometimes?
      //ASSERT(hg.connectivity(he) > 1, "Hyperedge is not a cut hyperedge!");
      for (const HypernodeID& pin : hg.pins(he)) {
        if (!visited[hg.originalNodeID(pin)]) {
          if (hg.partID(pin) == block_0) {
            if(scheduler.tryAquireNode(pin)){
              start_nodes_block_0.push_back(pin);
            }
          } else if (hg.partID(pin) == block_1) {
            if(scheduler.tryAquireNode(pin)){
              start_nodes_block_1.push_back(pin);
            }
          }
          visited.set(hg.originalNodeID(pin), true);
        }
      }
    }
    visited.reset();

    //TODO: partWeight of blocks can be influenced by other flow calculations. 
    // Use other metric or save and use them before a round?
    const HypernodeWeight max_part_weight_0 =
      std::max(((1.0 + std::min(alpha * context.partition.epsilon, 0.5))
                * context.partition.perfect_balance_part_weights[1]
                - context.partition.perfect_balance_part_weights[0]), 0.0); //hg.partWeight(block_1)
    const HypernodeWeight max_part_weight_1 =
      std::max(((1.0 + std::min(alpha * context.partition.epsilon, 0.5))
                * context.partition.perfect_balance_part_weights[0]
                - context.partition.perfect_balance_part_weights[0]), 0.0); //hg.partWeight(block_0)

    const HypernodeID num_nodes_block_0 = FlowRegionBuildPolicy<TypeTraits>::bfs(hg, flow_network,
                                                                     start_nodes_block_0,
                                                                     block_0,
                                                                     block_1,
                                                                     max_part_weight_0,
                                                                     visited,
                                                                     scheduler);
    if (num_nodes_block_0 == hg.partSize(block_0)) {
      // prevent blocks from becoming empty
      const HypernodeID last_hn_block_0 = hg.globalNodeID(*(flow_network.hypernodes().second - 1));
      flow_network.removeHypernode(hg, last_hn_block_0);
    }

    const HypernodeID num_nodes_block_1 = FlowRegionBuildPolicy<TypeTraits>::bfs(hg,
                                                                     flow_network,
                                                                     start_nodes_block_1,
                                                                     block_1,
                                                                     block_0,
                                                                     max_part_weight_1,
                                                                     visited,
                                                                     scheduler);
    if (num_nodes_block_1 == hg.partSize(block_1)) {
      // prevent blocks from becoming empty
      const HypernodeID last_hn_block_1 = hg.globalNodeID(*(flow_network.hypernodes().second - 1));
      flow_network.removeHypernode(hg, last_hn_block_1);
    }

    ASSERT([&]() {
        HypernodeWeight weight_block_0 = 0;
        HypernodeWeight weight_block_1 = 0;
        for (const HypernodeID& ogHn : flow_network.hypernodes()) {
          const HypernodeID& hn = hg.globalNodeID(ogHn);
          if (hg.partID(hn) == block_0) {
            weight_block_0 += hg.nodeWeight(hn);
          } else if (hg.partID(hn) == block_1) {
            weight_block_1 += hg.nodeWeight(hn);
          }
        }
        if (weight_block_0 > max_part_weight_0) {
          return false;
        }
        if (weight_block_1 > max_part_weight_1) {
          return false;
        }
        return true;
      } (), "Block size of one part is greater than the maximum allowed block size!");
  }
};
}  // namespace mt_kahypar
