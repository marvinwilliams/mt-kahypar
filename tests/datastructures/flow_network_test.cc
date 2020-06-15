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

#include "gmock/gmock.h"

#include "tests/datastructures/hypergraph_fixtures.h"
#include "mt-kahypar/datastructures/flow_network.h"
#include "mt-kahypar/definitions.h"
#include "mt-kahypar/partition/metrics.h"
#include "mt-kahypar/partition/refinement/flow/quotient_graph_block_scheduler.h"
#include "mt-kahypar/partition/refinement/flow/flow_region_build_policy.h"
#include "mt-kahypar/partition/refinement/flow/most_balanced_minimum_cut.h"

using ::testing::Test;

namespace mt_kahypar {
namespace ds {

using PartitionedHyperGraph = mt_kahypar::PartitionedHypergraph;

#define INCOMING(X) flowNetwork.mapToIncommingHyperedgeID(X)
#define OUTGOING(X) flowNetwork.mapToOutgoingHyperedgeID(X)
#define EDGE(X, C) std::make_pair(X, C)

typedef std::pair<NodeID, Capacity> edge;

using ds::FlowEdge;

template <typename FlowTypeTraits>
class AFlowNetworkTest : public ::testing::TestWithParam<std::pair<NodeID, std::set<edge> > >{
public:
    using Scheduler = typename FlowTypeTraits::Scheduler;
    using RegionBuildPolicy =  typename FlowTypeTraits::RegionBuildPolicy;
    using FlowNetwork = typename FlowTypeTraits::FlowNetwork;

  AFlowNetworkTest() :
    hg(HypergraphFactory::construct(TBBNumaArena::GLOBAL_TASK_GROUP,
      80 , 13, { { 0, 20},
                {1, 2, 22},
                {3, 8, 22},
                {20, 21, 26},
                {26, 31},
                {31, 32, 33},
                {1, 6, 7, 11, 12},
                {10, 11},
                {14, 50},
                {12, 13, 14, 17, 18},
                {16, 17},
                {4, 40, 41, 45},
                {50, 51, 56} })),
    hypergraph(),
    context(),
    visited(93),
    flow_network(),
    flow_network2(),
    scheduler() {
        context.partition.k = 4;
        context.partition.epsilon = 0.03;
        context.setupPartWeights(80);
        context.refinement.flow.alpha = 32;

        hypergraph = PartitionedHyperGraph(4, TBBNumaArena::GLOBAL_TASK_GROUP, hg);

        // Assign part ids
        for ( HypernodeID hn = 0; hn < 80; ++hn ) {
            hypergraph.setNodePart(hn, hn / 20);
        }
        hypergraph.initializeNumCutHyperedges();

        flow_network = std::make_unique<FlowNetwork>(80, 13, 106);
        flow_network2 = std::make_unique<FlowNetwork>(80, 13, 106);
        scheduler = std::make_unique<Scheduler>(hypergraph, context);
    }

    void buildFlowNetwork(PartitionedHyperGraph& hg,
                                    const Context& context,
                                    FlowNetwork& flow_network,
                                    std::vector<HyperedgeID>& cut_hes,
                                    const double alpha,
                                    const PartitionID block_0,
                                    const PartitionID block_1,
                                    kahypar::ds::FastResetFlagArray<>& visited,
                                    Scheduler & scheduler){
        RegionBuildPolicy::buildFlowNetwork(
                    hg, context, flow_network,
                    cut_hes, alpha, block_0, block_1,
                    visited, scheduler);
    }

  Hypergraph hg;
  PartitionedHyperGraph hypergraph;
  Context context;
  kahypar::ds::FastResetFlagArray<> visited;

  std::unique_ptr<FlowNetwork> flow_network;
  std::unique_ptr<FlowNetwork> flow_network2;
  std::unique_ptr<Scheduler> scheduler;
};

struct FlowTestOptTypeTraits{
    using Scheduler = OptScheduler;
    using RegionBuildPolicy = OptFlowRegionBuildPolicy;
    using FlowNetwork = ds::OptFlowNetwork<FlowTestOptTypeTraits>;
    using MostBalancedMinimumCut = OptMostBalancedMinimumCut<FlowTestOptTypeTraits>;
};

typedef ::testing::Types<FlowTestOptTypeTraits> OptConfig;

TYPED_TEST_CASE(AFlowNetworkTest, OptConfig);

TYPED_TEST(AFlowNetworkTest, FlowRegionTest) {
  this->scheduler->buildQuotientGraph();
  std::vector<HyperedgeID> cut_hes;
  for (const HyperedgeID& he : this->scheduler->blockPairCutHyperedges(0, 1)) {
      cut_hes.push_back(he);
  }
  this->flow_network->reset(0,1);
  this->buildFlowNetwork(
                  this->hypergraph, this->context, *this->flow_network,
                  cut_hes, this->context.refinement.flow.alpha, 0, 1,
                  this->visited, *this->scheduler);
  std::vector<HypernodeID> contained = {0,1, 2, 3, 6, 7, 8, 10, 11, 12, 20, 21, 22, 26, 31, 32, 33};
  std::vector<HypernodeID> not_contained = {4, 5, 9, 13, 14, 15, 16, 17, 18, 19, 23, 24, 25, 27, 28, 29, 30, 34, 35, 36, 37, 38, 39};
  for(HypernodeID node: contained){
    ASSERT_TRUE(this->flow_network->containsHypernode(node));
    ASSERT_TRUE(this->scheduler->isAquired(node));
  }
  for(HypernodeID node: not_contained){
    ASSERT_FALSE(this->flow_network->containsHypernode(node));
    ASSERT_FALSE(this->scheduler->isAquired(node));
  }

  std::vector<HyperedgeID> cut_hes2;
  for (const HyperedgeID& he : this->scheduler->blockPairCutHyperedges(0, 2)) {
      cut_hes2.push_back(he);
  }
  this->visited.reset();
  this->flow_network2->reset(0,2);
  this->buildFlowNetwork(
                  this->hypergraph, this->context, *this->flow_network2,
                  cut_hes2, this->context.refinement.flow.alpha, 0, 2,
                  this->visited, *this->scheduler);
  contained = {4, 13, 14, 16, 17, 18, 40, 41, 45, 50, 51, 56};
  not_contained = {5, 9, 15, 19};
  for(HypernodeID node: contained){
    ASSERT_TRUE(this->flow_network2->containsHypernode(node));
    ASSERT_TRUE(this->scheduler->isAquired(node));
  }
  for(HypernodeID node: not_contained){
    ASSERT_FALSE(this->flow_network2->containsHypernode(node));
    ASSERT_FALSE(this->scheduler->isAquired(node));
  }
}

TYPED_TEST(AFlowNetworkTest, Release) {
  this->scheduler->buildQuotientGraph();
  std::vector<HyperedgeID> cut_hes;
  for (const HyperedgeID& he : this->scheduler->blockPairCutHyperedges(0, 1)) {
      cut_hes.push_back(he);
  }
  this->flow_network->reset(0,1);
  this->buildFlowNetwork(
                  this->hypergraph, this->context, *this->flow_network,
                  cut_hes, this->context.refinement.flow.alpha, 0, 1,
                  this->visited, *this->scheduler);
  std::vector<HyperedgeID> cut_hes2;
  for (const HyperedgeID& he : this->scheduler->blockPairCutHyperedges(0, 2)) {
      cut_hes2.push_back(he);
  }
  this->visited.reset();
  this->flow_network2->reset(0,2);
  this->buildFlowNetwork(
                  this->hypergraph, this->context, *this->flow_network2,
                  cut_hes2, this->context.refinement.flow.alpha, 0, 2,
                  this->visited, *this->scheduler);

  std::vector<HypernodeID> aquired = {0,1, 2, 3, 6, 7, 8, 10, 11, 12, 20, 21, 22, 26, 31, 32, 33,
                                         4, 13, 14, 16, 17, 18, 40, 41, 45, 50, 51, 56};
  for(HypernodeID node: aquired){
    ASSERT_TRUE(this->scheduler->isAquired(node));
  }
  this->flow_network->release(this->hypergraph, 0, 1, *this->scheduler);
  this->flow_network2->release(this->hypergraph, 0, 2, *this->scheduler);
  for(HypernodeID node: aquired){
    ASSERT_FALSE(this->scheduler->isAquired(node));
  }
  std::pair<HypernodeWeight, HypernodeWeight> weights = this->scheduler->get_aquired_part_weight(0, 1);
  ASSERT_EQ(weights.first, 0);
  ASSERT_EQ(weights.second, 0);

  weights = this->scheduler->get_aquired_part_weight(0, 2);
  ASSERT_EQ(weights.first, 0);
  ASSERT_EQ(weights.second, 0);
}

TYPED_TEST(AFlowNetworkTest, FlowNetWorkBuild) {
  this->scheduler->buildQuotientGraph();
  std::vector<HyperedgeID> cut_hes;
  for (const HyperedgeID& he : this->scheduler->blockPairCutHyperedges(0, 1)) {
      cut_hes.push_back(he);
  }
  this->flow_network->reset(0,1);
  this->buildFlowNetwork(
                  this->hypergraph, this->context, *this->flow_network,
                  cut_hes, this->context.refinement.flow.alpha, 0, 1,
                  this->visited, *this->scheduler);

  std::vector<HyperedgeID> cut_hes2;
  for (const HyperedgeID& he : this->scheduler->blockPairCutHyperedges(0, 2)) {
      cut_hes2.push_back(he);
  }
  this->visited.reset();
  this->flow_network2->reset(0,2);
  this->buildFlowNetwork(
                  this->hypergraph, this->context, *this->flow_network2,
                  cut_hes2, this->context.refinement.flow.alpha, 0, 2,
                  this->visited, *this->scheduler);

  this->flow_network->build(this->hypergraph, this->context, 0, 1, *this->scheduler);

  ASSERT_EQ(this->flow_network->numNodes(), 18);
  ASSERT_EQ(this->flow_network->numEdges(), 26);
  ASSERT_EQ(this->flow_network->numUndirectedEdges(), 3);
  ASSERT_EQ(this->flow_network->totalWeightHyperedges(), 8);
  ASSERT_TRUE(this->flow_network->isRemovedHypernode(21));

  this->flow_network2->build(this->hypergraph, this->context, 0, 2, *this->scheduler);

  ASSERT_EQ(this->flow_network2->numNodes(), 11);
  ASSERT_EQ(this->flow_network2->numEdges(), 12);
  ASSERT_EQ(this->flow_network2->numUndirectedEdges(), 2);
  ASSERT_EQ(this->flow_network2->totalWeightHyperedges(), 4);
  ASSERT_TRUE(this->flow_network2->isRemovedHypernode(51));
}

TYPED_TEST(AFlowNetworkTest, SourcesAndSinks) {
  this->scheduler->buildQuotientGraph();
  std::vector<HyperedgeID> cut_hes;
  for (const HyperedgeID& he : this->scheduler->blockPairCutHyperedges(0, 1)) {
      cut_hes.push_back(he);
  }
  this->flow_network->reset(0,1);
  this->buildFlowNetwork(
                  this->hypergraph, this->context, *this->flow_network,
                  cut_hes, this->context.refinement.flow.alpha, 0, 1,
                  this->visited, *this->scheduler);

  std::vector<HyperedgeID> cut_hes2;
  for (const HyperedgeID& he : this->scheduler->blockPairCutHyperedges(0, 2)) {
      cut_hes2.push_back(he);
  }
  this->visited.reset();
  this->flow_network2->reset(0,2);
  this->buildFlowNetwork(
                  this->hypergraph, this->context, *this->flow_network2,
                  cut_hes2, this->context.refinement.flow.alpha, 0, 2,
                  this->visited, *this->scheduler);

  this->flow_network->build(this->hypergraph, this->context, 0, 1, *this->scheduler);
  this->flow_network2->build(this->hypergraph, this->context, 0, 2, *this->scheduler);

  ASSERT_EQ(this->flow_network->numSources(), 1);
  ASSERT_TRUE(this->flow_network->isIdSource(89));

  ASSERT_EQ(this->flow_network2->numSources(), 1);
  ASSERT_TRUE(this->flow_network2->isIdSource(89));
}

TYPED_TEST(AFlowNetworkTest, IncidentEdges) {
  this->scheduler->buildQuotientGraph();
  std::vector<HyperedgeID> cut_hes;
  for (const HyperedgeID& he : this->scheduler->blockPairCutHyperedges(0, 1)) {
      cut_hes.push_back(he);
  }
  this->flow_network->reset(0,1);
  this->buildFlowNetwork(
                  this->hypergraph, this->context, *this->flow_network,
                  cut_hes, this->context.refinement.flow.alpha, 0, 1,
                  this->visited, *this->scheduler);

  std::vector<HyperedgeID> cut_hes2;
  for (const HyperedgeID& he : this->scheduler->blockPairCutHyperedges(0, 2)) {
      cut_hes2.push_back(he);
  }
  this->visited.reset();
  this->flow_network2->reset(0,2);
  this->buildFlowNetwork(
                  this->hypergraph, this->context, *this->flow_network2,
                  cut_hes2, this->context.refinement.flow.alpha, 0, 2,
                  this->visited, *this->scheduler);

  this->flow_network->build(this->hypergraph, this->context, 0, 1, *this->scheduler);
  this->flow_network2->build(this->hypergraph, this->context, 0, 2, *this->scheduler);

  std::set<size_t> targets = {86, 99, 89};
  for(auto e:this->flow_network->incidentEdges(12)){
    ASSERT_TRUE(targets.find(e.target) != targets.end());
  }
  targets = {50, 89};
  for(auto e:this->flow_network2->incidentEdges(14)){
    ASSERT_TRUE(targets.find(e.target) != targets.end());
  }
}

TYPED_TEST(AFlowNetworkTest, AquiredPartWeight) {
  this->scheduler->buildQuotientGraph();
  std::vector<HyperedgeID> cut_hes;
  for (const HyperedgeID& he : this->scheduler->blockPairCutHyperedges(0, 1)) {
      cut_hes.push_back(he);
  }
  this->flow_network->reset(0,1);
  this->buildFlowNetwork(
                  this->hypergraph, this->context, *this->flow_network,
                  cut_hes, this->context.refinement.flow.alpha, 0, 1,
                  this->visited, *this->scheduler);

  std::vector<HyperedgeID> cut_hes2;
  for (const HyperedgeID& he : this->scheduler->blockPairCutHyperedges(0, 2)) {
      cut_hes2.push_back(he);
  }
  this->visited.reset();
  this->flow_network2->reset(0,2);
  this->buildFlowNetwork(
                  this->hypergraph, this->context, *this->flow_network2,
                  cut_hes2, this->context.refinement.flow.alpha, 0, 2,
                  this->visited, *this->scheduler);

  this->flow_network->build(this->hypergraph, this->context, 0, 1, *this->scheduler);
  this->flow_network2->build(this->hypergraph, this->context, 0, 2, *this->scheduler);

  auto partweights = this->scheduler->get_aquired_part_weight(0, 1);
  ASSERT_EQ(partweights.first, 10);
  ASSERT_EQ(partweights.second, 7);

  partweights = this->scheduler->get_aquired_part_weight(0, 2);
  ASSERT_EQ(partweights.first, 6);
  ASSERT_EQ(partweights.second, 6);
}

}  //namespace ds
}  // namespace mt_kahypar
