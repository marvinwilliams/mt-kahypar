/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2020 Tobias Heuer <tobias.heuer@kit.edu>
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
#include <atomic>

#include "mt-kahypar/datastructures/dynamic_adjacency_array.h"
#include "mt-kahypar/parallel/atomic_wrapper.h"

using ::testing::Test;

namespace mt_kahypar {
namespace ds {

void verifyNeighbors(const HypernodeID u,
                     const HypernodeID num_nodes,
                     const DynamicAdjacencyArray& incident_edges,
                     const std::set<HypernodeID>& _expected_neighbors,
                     bool strict = false) {
  size_t num_neighbors = 0;
  size_t degree = 0;
  std::vector<bool> actual_neighbors(num_nodes, false);
  for ( const HyperedgeID& he : incident_edges.incidentEdges(u) ) {
    const HypernodeID neighbor = incident_edges.edge(he).target;
    ASSERT_NE(_expected_neighbors.find(neighbor), _expected_neighbors.end())
      << "Vertex " << neighbor << " should not be neighbor of vertex " << u;
    ASSERT_EQ(u, incident_edges.edge(he).source)
      << "Source of " << he << " (target: " << incident_edges.edge(he).target << ") should be "
      << u << " but is " << incident_edges.edge(he).source;
    ASSERT_TRUE(!strict || !actual_neighbors[neighbor])
      << "Vertex " << u << " contain duplicate edge with target " << neighbor;
    if (!actual_neighbors[neighbor]) {
      ++num_neighbors;
    }
    ++degree;
    actual_neighbors[neighbor] = true;
  }
  ASSERT_EQ(num_neighbors, _expected_neighbors.size());
  ASSERT_EQ(degree, incident_edges.nodeDegree(u));
  ASSERT_TRUE(!strict || num_neighbors == degree);
}

kahypar::ds::FastResetFlagArray<> createFlagArray(const HypernodeID num_nodes,
                                                  const std::vector<HypernodeID>& contained_nodes) {
  kahypar::ds::FastResetFlagArray<> flag_array(num_nodes);
  for ( const HypernodeID& node : contained_nodes ) {
    flag_array.set(node, true);
  }
  return flag_array;
}

TEST(ADynamicAdjacencyArray, VerifyInitialNeighborsOfEachVertex) {
  DynamicAdjacencyArray incident_edges(
    7, {{1, 2}, {2, 3}, {1, 4}, {4, 5}, {4, 6}, {5, 6}});
  verifyNeighbors(0, 7, incident_edges, { });
  verifyNeighbors(1, 7, incident_edges, { 2, 4 });
  verifyNeighbors(2, 7, incident_edges, { 1, 3 });
  verifyNeighbors(3, 7, incident_edges, { 2 });
  verifyNeighbors(4, 7, incident_edges, { 1, 5, 6 });
  verifyNeighbors(5, 7, incident_edges, { 4, 6 });
  verifyNeighbors(6, 7, incident_edges, { 4, 5 });
}

TEST(ADynamicAdjacencyArray, ContractTwoVertices1) {
  DynamicAdjacencyArray incident_edges(
    7, {{1, 2}, {2, 3}, {1, 4}, {4, 5}, {4, 6}, {5, 6}});
  auto flag_array = createFlagArray(7, { });
  incident_edges.contract(0, 1, flag_array);
  verifyNeighbors(0, 7, incident_edges, { 2, 4 });
}

TEST(ADynamicAdjacencyArray, ContractTwoVertices2) {
  DynamicAdjacencyArray incident_edges(
    7, {{1, 2}, {2, 3}, {1, 4}, {4, 5}, {4, 6}, {5, 6}});
  auto flag_array = createFlagArray(7, { });
  incident_edges.contract(1, 2, flag_array);
  verifyNeighbors(1, 7, incident_edges, { 3, 4 });
}

TEST(ADynamicAdjacencyArray, ContractTwoVertices3) {
  DynamicAdjacencyArray incident_edges(
    7, {{1, 2}, {2, 3}, {1, 4}, {4, 5}, {4, 6}, {5, 6}});
  auto flag_array = createFlagArray(7, { });
  incident_edges.contract(1, 5, flag_array);
  verifyNeighbors(1, 7, incident_edges, { 2, 4, 6 });
  verifyNeighbors(2, 7, incident_edges, { 1, 3 });
  verifyNeighbors(4, 7, incident_edges, { 1, 6 });
  verifyNeighbors(6, 7, incident_edges, { 1, 4 });
}

TEST(ADynamicAdjacencyArray, ContractSeveralVertices1) {
  DynamicAdjacencyArray incident_edges(
    7, {{1, 2}, {2, 3}, {1, 4}, {4, 5}, {4, 6}, {5, 6}});
  auto flag_array = createFlagArray(7, { });
  incident_edges.contract(0, 1, flag_array);
  incident_edges.contract(0, 2, flag_array);
  verifyNeighbors(0, 7, incident_edges, { 3, 4 });
}

TEST(ADynamicAdjacencyArray, ContractSeveralVertices2) {
  DynamicAdjacencyArray incident_edges(
    7, {{1, 2}, {2, 3}, {1, 4}, {4, 5}, {4, 6}, {5, 6}});
  auto flag_array = createFlagArray(7, { });
  incident_edges.contract(1, 0, flag_array);
  incident_edges.contract(1, 2, flag_array);
  verifyNeighbors(1, 7, incident_edges, { 3, 4 });
  incident_edges.contract(4, 5, flag_array);
  incident_edges.contract(4, 6, flag_array);
  verifyNeighbors(4, 7, incident_edges, { 1 });
  incident_edges.contract(1, 3, flag_array);
  incident_edges.contract(1, 4, flag_array);
  verifyNeighbors(1, 7, incident_edges, { });
}

TEST(ADynamicAdjacencyArray, UncontractTwoVertices1) {
  DynamicAdjacencyArray incident_edges(
    7, {{1, 2}, {2, 3}, {1, 4}, {4, 5}, {4, 6}, {5, 6}});
  auto flag_array = createFlagArray(7, { });
  incident_edges.contract(1, 2, flag_array);
  incident_edges.uncontract(1, 2);
  verifyNeighbors(1, 7, incident_edges, { 2, 4 });
  verifyNeighbors(2, 7, incident_edges, { 1, 3 });
}

TEST(ADynamicAdjacencyArray, UncontractTwoVertices2) {
  DynamicAdjacencyArray incident_edges(
    7, {{1, 2}, {2, 3}, {1, 4}, {4, 5}, {4, 6}, {5, 6}});
  auto flag_array = createFlagArray(7, { });
  incident_edges.contract(1, 5, flag_array);
  incident_edges.uncontract(1, 5);
  verifyNeighbors(1, 7, incident_edges, { 2, 4 });
  verifyNeighbors(5, 7, incident_edges, { 4, 6 });
  verifyNeighbors(2, 7, incident_edges, { 1, 3 });
  verifyNeighbors(4, 7, incident_edges, { 1, 5, 6 });
  verifyNeighbors(6, 7, incident_edges, { 4, 5 });
}

TEST(ADynamicAdjacencyArray, UncontractSeveralVertices1) {
  DynamicAdjacencyArray incident_edges(
    7, {{1, 2}, {2, 3}, {1, 4}, {4, 5}, {4, 6}, {5, 6}});
  auto flag_array = createFlagArray(7, { });
  incident_edges.contract(1, 2, flag_array);
  incident_edges.contract(1, 0, flag_array);
  incident_edges.contract(4, 5, flag_array);
  incident_edges.contract(4, 6, flag_array);
  incident_edges.contract(4, 3, flag_array);
  incident_edges.contract(4, 1, flag_array);
  verifyNeighbors(4, 7, incident_edges, { });
  incident_edges.uncontract(4, 1);
  verifyNeighbors(1, 7, incident_edges, { 4 });
  verifyNeighbors(4, 7, incident_edges, { 1 });
  incident_edges.uncontract(4, 3);
  verifyNeighbors(4, 7, incident_edges, { 1 });
  verifyNeighbors(3, 7, incident_edges, { 1 });
  incident_edges.uncontract(4, 6);
  verifyNeighbors(4, 7, incident_edges, { 1, 6 });
  verifyNeighbors(6, 7, incident_edges, { 4 });
  incident_edges.uncontract(4, 5);
  verifyNeighbors(4, 7, incident_edges, { 1, 5, 6 });
  verifyNeighbors(5, 7, incident_edges, { 4, 6 });
  incident_edges.uncontract(1, 0);
  verifyNeighbors(0, 7, incident_edges, { });
  verifyNeighbors(1, 7, incident_edges, { 3, 4 });
  incident_edges.uncontract(1, 2);
  verifyNeighbors(0, 7, incident_edges, { });
  verifyNeighbors(1, 7, incident_edges, { 2, 4 });
  verifyNeighbors(2, 7, incident_edges, { 1, 3 });
  verifyNeighbors(3, 7, incident_edges, { 2 });
  verifyNeighbors(4, 7, incident_edges, { 1, 5, 6 });
  verifyNeighbors(5, 7, incident_edges, { 4, 6 });
  verifyNeighbors(6, 7, incident_edges, { 4, 5 });
}

TEST(ADynamicAdjacencyArray, UncontractSeveralVertices2) {
  DynamicAdjacencyArray incident_edges(
    7, {{1, 2}, {2, 3}, {1, 4}, {4, 5}, {4, 6}, {5, 6}});
  auto flag_array = createFlagArray(7, { });
  incident_edges.contract(3, 1, flag_array);
  incident_edges.contract(3, 4, flag_array);
  incident_edges.contract(5, 6, flag_array);
  incident_edges.contract(3, 5, flag_array);
  incident_edges.contract(0, 2, flag_array);
  incident_edges.contract(0, 3, flag_array);
  verifyNeighbors(0, 7, incident_edges, { });
  incident_edges.uncontract(0, 3);
  verifyNeighbors(0, 7, incident_edges, { 3 });
  verifyNeighbors(3, 7, incident_edges, { 0 });
  incident_edges.uncontract(0, 2);
  verifyNeighbors(0, 7, incident_edges, { });
  verifyNeighbors(2, 7, incident_edges, { 3 });
  incident_edges.uncontract(3, 5);
  verifyNeighbors(3, 7, incident_edges, { 2, 5 });
  verifyNeighbors(5, 7, incident_edges, { 3 });
  incident_edges.uncontract(5, 6);
  verifyNeighbors(5, 7, incident_edges, { 3, 6 });
  verifyNeighbors(6, 7, incident_edges, { 3, 5 });
  incident_edges.uncontract(3, 4);
  verifyNeighbors(3, 7, incident_edges, { 2, 4 });
  verifyNeighbors(4, 7, incident_edges, { 3, 5, 6 });
  incident_edges.uncontract(3, 1);
  verifyNeighbors(0, 7, incident_edges, { });
  verifyNeighbors(1, 7, incident_edges, { 2, 4 });
  verifyNeighbors(2, 7, incident_edges, { 1, 3 });
  verifyNeighbors(3, 7, incident_edges, { 2 });
  verifyNeighbors(4, 7, incident_edges, { 1, 5, 6 });
  verifyNeighbors(5, 7, incident_edges, { 4, 6 });
  verifyNeighbors(6, 7, incident_edges, { 4, 5 });
}

TEST(ADynamicAdjacencyArray, RemovesParrallelEdges1) {
  DynamicAdjacencyArray incident_edges(
    7, {{1, 2}, {2, 3}, {1, 4}, {4, 5}, {4, 6}, {5, 6}});
  auto flag_array = createFlagArray(7, { });
  incident_edges.contract(2, 4, flag_array);
  incident_edges.removeParallelEdges();
  verifyNeighbors(1, 7, incident_edges, { 2 }, true);
  verifyNeighbors(2, 7, incident_edges, { 1, 3, 5, 6 }, true);
}


// using OwnershipVector = parallel::scalable_vector<SpinLock>;

// template<typename F, typename K>
// void executeParallel(const F& f1, const K& f2) {
//   std::atomic<size_t> cnt(0);
//   tbb::parallel_invoke([&] {
//     ++cnt;
//     while ( cnt < 2 ) { }
//     f1();
//   }, [&] {
//     ++cnt;
//     while ( cnt < 2 ) { }
//     f2();
//   });
// }


}  // namespace ds
}  // namespace mt_kahypar
