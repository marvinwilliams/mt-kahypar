/*******************************************************************************
 * This file is part of KaHyPar.
 *
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

#include "gmock/gmock.h"

#include "tests/datastructures/hypergraph_fixtures.h"
#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/partition/refinement/flow/quotient_graph_block_scheduler.h"


namespace mt_kahypar {

using ::testing::Test;

template<typename Scheduler>
class AOptScheduler : public Test {
 public:
  AOptScheduler() :
    hg(HypergraphFactory::construct(TBBNumaArena::GLOBAL_TASK_GROUP,
      32 , 8, { { 0, 1,  2,  3,  8,  9,  10,  11 },
                {0, 20},
                {0, 24},
                { 4, 5, 6, 7, 12, 13, 14, 15 },
                {4, 28},
                {12, 28},
                {16, 20},
                {24,28} })),
    hypergraph(),
    context(),
    scheduler() {

    context.partition.k = 8;
    context.partition.epsilon = 0.0;
    context.setupPartWeights(32);

    hypergraph = PartitionedHypergraph<>(8, TBBNumaArena::GLOBAL_TASK_GROUP, hg);

    // Assign part ids
    for ( HypernodeID hn = 0; hn < 32; ++hn ) {
      hypergraph.setNodePart(hn, hn / 4);
    }
    hypergraph.initializeNumCutHyperedges();

    scheduler = std::make_unique<Scheduler>(hypergraph, context);
  }

  Hypergraph hg;
  PartitionedHypergraph<> hypergraph;
  Context context;

  std::unique_ptr<Scheduler> scheduler;
};

typedef ::testing::Types<OptScheduler> OptConfig;

TYPED_TEST_CASE(AOptScheduler, OptConfig);


TYPED_TEST(AOptScheduler, GetInitialEdgesOpt) {
  this->scheduler->buildQuotientGraph();
  std::vector<edge> initial_edges = this->scheduler->getInitialParallelEdges();
  std::vector<edge> right_init_edges = { std::make_pair(0, 2), std::make_pair(6, 7),
                                         std::make_pair(4, 5), std::make_pair(1, 3),
                                         std::make_pair(3, 7), std::make_pair(0, 5),
                                         std::make_pair(1, 7),std::make_pair(0, 6)} ;
  for (size_t i = 0; i < initial_edges.size(); i++){
    ASSERT_EQ(initial_edges[i].first, right_init_edges[i].first);
    ASSERT_EQ(initial_edges[i].second, right_init_edges[i].second);
  }
}

TYPED_TEST(AOptScheduler, HasCorrectCutHyperedges) {
  this->scheduler->buildQuotientGraph();

  for (const auto& e : this->scheduler->blockPairCutHyperedges(0, 2)) {
    ASSERT_EQ(e, 0);
  }
  for (const auto& e : this->scheduler->blockPairCutHyperedges(0, 5)) {
    ASSERT_EQ(e, 1);
  }
  for (const auto& e : this->scheduler->blockPairCutHyperedges(0, 6)) {
    ASSERT_EQ(e, 2);
  }
  for (const auto& e : this->scheduler->blockPairCutHyperedges(1, 3)) {
    ASSERT_EQ(e, 3);
  }
  for (const auto& e : this->scheduler->blockPairCutHyperedges(1, 7)) {
    ASSERT_EQ(e, 4);
  }
  for (const auto& e : this->scheduler->blockPairCutHyperedges(3, 7)) {
    ASSERT_EQ(e, 5);
  }
  for (const auto& e : this->scheduler->blockPairCutHyperedges(4, 5)) {
    ASSERT_EQ(e, 6);
  }
  for (const auto& e : this->scheduler->blockPairCutHyperedges(6, 7)) {
    ASSERT_EQ(e, 7);
  }
}

TYPED_TEST(AOptScheduler, HasCorrectCutHyperedgesAfterMove) {
  this->scheduler->buildQuotientGraph();

  this->scheduler->changeNodePart(28, 7, 4);

  for (const auto& e : this->scheduler->blockPairCutHyperedges(0, 2)) {
    ASSERT_EQ(e, 0);
  }
  for (const auto& e : this->scheduler->blockPairCutHyperedges(0, 5)) {
    ASSERT_EQ(e, 1);
  }
  for (const auto& e : this->scheduler->blockPairCutHyperedges(0, 6)) {
    ASSERT_EQ(e, 2);
  }
  for (const auto& e : this->scheduler->blockPairCutHyperedges(1, 3)) {
    ASSERT_EQ(e, 3);
  }
  size_t idx = 0;
  for (const auto& e : this->scheduler->blockPairCutHyperedges(1, 7)) {
    unused(e);
    idx++;
  }
  ASSERT_EQ(idx, 0);
  idx = 0;
  for (const auto& e : this->scheduler->blockPairCutHyperedges(3, 7)) {
    unused(e);
    idx ++;
  }
  ASSERT_EQ(idx, 0);
  for (const auto& e : this->scheduler->blockPairCutHyperedges(4, 5)) {
    ASSERT_EQ(e, 6);
  }
  idx = 0;
  for (const auto& e : this->scheduler->blockPairCutHyperedges(6, 7)) {
    unused(e);
    idx++;
  }
  ASSERT_EQ(idx, 0);

  for (const auto& e : this->scheduler->blockPairCutHyperedges(1, 4)) {
    ASSERT_EQ(e, 4);
  }
  for (const auto& e : this->scheduler->blockPairCutHyperedges(3, 4)) {
    ASSERT_EQ(e, 5);
  }
  for (const auto& e : this->scheduler->blockPairCutHyperedges(4, 6)) {
    ASSERT_EQ(e, 7);
  }


}

TYPED_TEST(AOptScheduler, AquireNodeTest) {
  size_t block_idx = 1;
  HypernodeID node = 0;
  bool success = this->scheduler->tryAquireNode(node, block_idx);
  ASSERT_TRUE(success);

  success = this->scheduler->tryAquireNode(node, block_idx);
  ASSERT_FALSE(success);

  bool aquired = this->scheduler->isAquired(node);
  ASSERT_TRUE(aquired);

  this->scheduler->releaseNode(node);
  aquired = this->scheduler->isAquired(node);
  ASSERT_FALSE(aquired);
}

TYPED_TEST(AOptScheduler, BlockOverlapTest) {
  size_t block_idx = 1;
  HypernodeID node = 0;
  this->scheduler->tryAquireNode(node, block_idx);
  bool overlap = this->scheduler->is_block_overlap(node, 0 , 2);
  ASSERT_TRUE(overlap);

  overlap = this->scheduler->is_block_overlap(node, 1 , 2);
  ASSERT_TRUE(overlap);

  overlap = this->scheduler->is_block_overlap(node, 2 , 3);
  ASSERT_FALSE(overlap);

  this->scheduler->releaseNode(node);
  block_idx = 3 * this->context.partition.k + 4;
  this->scheduler->tryAquireNode(node, block_idx);
  overlap = this->scheduler->is_block_overlap(node, 2 , 3);
  ASSERT_TRUE(overlap);
}

TYPED_TEST(AOptScheduler, ActiveBlocks) {
  this->scheduler->buildQuotientGraph();
  ASSERT_EQ(this->context.partition.k , this->scheduler->getNumberOfActiveBlocks());

  std::vector<edge> initial_edges = this->scheduler->getInitialParallelEdges();
  ASSERT_EQ(0 , this->scheduler->getNumberOfActiveBlocks());

  this->scheduler->setActiveBlock(0, true);
  this->scheduler->setActiveBlock(1, true);
  ASSERT_EQ(2 , this->scheduler->getNumberOfActiveBlocks());
  this->scheduler->setActiveBlock(0, false);
  ASSERT_EQ(1 , this->scheduler->getNumberOfActiveBlocks());
}

TYPED_TEST(AOptScheduler, InitBlockWeight) {
  this->scheduler->init_block_weights();
  for (int i = 0; i < this->context.partition.k; i++){
    ASSERT_EQ(4 ,this->scheduler->get_not_aquired_weight(i, (i + 1)% this->context.partition.k));
  }
}

TYPED_TEST(AOptScheduler, AquireBlockWeight) {
  this->scheduler->init_block_weights();
  size_t block = 0;
  size_t other = 1;
  size_t amount = 3;
  this->scheduler->aquire_block_weight(block, other, amount);
  ASSERT_EQ(1 ,this->scheduler->get_not_aquired_weight(block, other));

  this->scheduler->release_block_weight(block, other, 2);
  ASSERT_EQ(3 ,this->scheduler->get_not_aquired_weight(block, other));
  other = 2;
  ASSERT_EQ(3 ,this->scheduler->get_not_aquired_weight(block, other));
}

template<typename Scheduler>
class AMatchScheduler : public Test {
 public:
  AMatchScheduler() :
    hg(HypergraphFactory::construct(TBBNumaArena::GLOBAL_TASK_GROUP,
      32 , 8, { { 0, 1,  2,  3,  8,  9,  10,  11 },
                {0, 20},
                {0, 24},
                { 4, 5, 6, 7, 12, 13, 14, 15 },
                {4, 28},
                {12, 28},
                {16, 20},
                {24,28} })),
    hypergraph(),
    context(),
    scheduler() {

    context.partition.k = 8;
    context.partition.epsilon = 0.0;
    context.setupPartWeights(32);

    hypergraph = PartitionedHypergraph<>(8, TBBNumaArena::GLOBAL_TASK_GROUP, hg);

    // Assign part ids
    for ( HypernodeID hn = 0; hn < 32; ++hn ) {
      hypergraph.setNodePart(hn, hn / 4);
    }
    hypergraph.initializeNumCutHyperedges();

    scheduler = std::make_unique<Scheduler>(hypergraph, context);
  }

  Hypergraph hg;
  PartitionedHypergraph<> hypergraph;
  Context context;

  std::unique_ptr<Scheduler> scheduler;
};

typedef ::testing::Types<MatchingScheduler> MatchConfig;

TYPED_TEST_CASE(AMatchScheduler, MatchConfig);


TYPED_TEST(AMatchScheduler, GetInitialEdgesMatch) {
  this->scheduler->buildQuotientGraph();
  std::vector<edge> initial_edges = this->scheduler->getInitialParallelEdges();
  std::vector<edge> right_init_edges = { std::make_pair(0, 2), std::make_pair(6, 7) ,
                                         std::make_pair(4, 5), std::make_pair(1, 3)};
  for (size_t i = 0; i < initial_edges.size(); i++){
    ASSERT_EQ(initial_edges[i].first, right_init_edges[i].first);
    ASSERT_EQ(initial_edges[i].second, right_init_edges[i].second);
  }
}


}//namespace Mt-KaHyPar