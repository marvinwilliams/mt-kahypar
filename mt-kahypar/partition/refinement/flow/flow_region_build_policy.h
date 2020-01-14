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
    using HyperGraph = typename TypeTraits::HyperGraph;


 public:
  template <class Network = Mandatory>
  static inline HypernodeID bfs(HyperGraph& hg,
                                Network& flow_network,
                                std::vector<HypernodeID>& start_nodes,
                                const PartitionID part,
                                const HypernodeWeight max_part_weight,
                                kahypar::ds::FastResetFlagArray<>& visited) {
    visited.reset();
    utils::Randomize::instance().shuffleVector(start_nodes);
    std::queue<HypernodeID> Q;
    HypernodeWeight queue_weight = 0;
    for (const HypernodeID& hn : start_nodes) {
      if (queue_weight + hg.nodeWeight(hn) <= max_part_weight) {
        Q.push(hn);
        queue_weight += hg.nodeWeight(hn);
        // TODO(reister): Note, in case we have several numa nodes, than node ids
        // are not consecutive. The first 16 bit are reserved for the numa node id
        // and the remaining 48 bits for the vertex id. However, you can map those
        // ids to a consecutive range by calling hg.originalNodeID(hn). Your line
        // would than look like this:
        // visited.set(hg.originalNodeID(hn), true);
        // You can test your code locally on a numa architecture by enabling
        // USE_HARDWARE_MOCK in definitions.h. Please, reconsider your code on all
        // places where you access an array with a vertex id and make sure that if
        // USE_HARDWARE_MOCK is enabled no segmentation fault occurs.
        visited.set(hg.originalNodeID(hn), true);
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
              Q.push(pin);
              queue_weight += hg.nodeWeight(pin);
              visited.set(hg.originalNodeID(pin), true);
            }
          }
          visited.set(num_hypernodes + hg.originalEdgeID(he), true);
        }
      }
    }
    return num_hypernodes_added;
  }
};

template< typename TypeTraits>
class CutBuildPolicy : public FlowRegionBuildPolicy<TypeTraits> {
  private:
    using HyperGraph = typename TypeTraits::HyperGraph;
 public:
  template <class Network = Mandatory>
  inline static void buildFlowNetwork(HyperGraph& hg,
                                      const Context& context,
                                      Network& flow_network,
                                      std::vector<HyperedgeID>& cut_hes,
                                      const double alpha,
                                      const PartitionID block_0,
                                      const PartitionID block_1,
                                      kahypar::ds::FastResetFlagArray<>& visited) {
    visited.reset();
    std::vector<HypernodeID> start_nodes_block_0;
    std::vector<HypernodeID> start_nodes_block_1;
    for (const HyperedgeID he : cut_hes) {
      ASSERT(hg.connectivity(he) > 1, "Hyperedge is not a cut hyperedge!");
      for (const HypernodeID& pin : hg.pins(he)) {
        if (!visited[hg.originalNodeID(pin)]) {
          if (hg.partID(pin) == block_0) {
            start_nodes_block_0.push_back(pin);
          } else if (hg.partID(pin) == block_1) {
            start_nodes_block_1.push_back(pin);
          }
          visited.set(hg.originalNodeID(pin), true);
        }
      }
    }
    visited.reset();

    const HypernodeWeight max_part_weight_0 =
      std::max(((1.0 + std::min(alpha * context.partition.epsilon, 0.5))
                * context.partition.perfect_balance_part_weights[1]
                - hg.localPartWeight(block_1)), 0.0);
    const HypernodeWeight max_part_weight_1 =
      std::max(((1.0 + std::min(alpha * context.partition.epsilon, 0.5))
                * context.partition.perfect_balance_part_weights[0]
                - hg.localPartWeight(block_0)), 0.0);

    const HypernodeID num_nodes_block_0 = FlowRegionBuildPolicy<TypeTraits>::bfs(hg, flow_network,
                                                                     start_nodes_block_0,
                                                                     block_0,
                                                                     max_part_weight_0,
                                                                     visited);
    if (num_nodes_block_0 == hg.localPartSize(block_0)) {
      // prevent blocks from becoming empty
      const HypernodeID last_hn_block_0 = hg.globalNodeID(*(flow_network.hypernodes().second - 1));
      flow_network.removeHypernode(hg, last_hn_block_0);
    }

    const HypernodeID num_nodes_block_1 = FlowRegionBuildPolicy<TypeTraits>::bfs(hg,
                                                                     flow_network,
                                                                     start_nodes_block_1,
                                                                     block_1,
                                                                     max_part_weight_1,
                                                                     visited);
    if (num_nodes_block_1 == hg.localPartSize(block_1)) {
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
