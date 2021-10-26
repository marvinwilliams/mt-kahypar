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

#include <atomic>
#include <cstdlib>
#include <vector>
#include <mt-kahypar/macros.h>

#include "gmock/gmock.h"

#include "mt-kahypar/utils/randomize.h"
#include "mt-kahypar/datastructures/perfect_hash_table.h"

using ::testing::Test;

#define P(X,Y) std::make_pair(X,Y)

namespace mt_kahypar {
namespace ds {

template<typename K, typename V>
void construct(PerfectHashTable<K, V>& pht, const std::vector<std::pair<K,V>>& elements) {
  pht.clear();
  for ( const auto& element : elements ) {
    pht.add(element.first, element.second);
  }
  pht.construct();
}

template<typename K, typename V>
PerfectHashTable<K, V> construct(const std::vector<std::pair<K,V>>& elements, const K universe_size) {
  PerfectHashTable<K,V> pht(universe_size);
  construct(pht, elements);
  return pht;
}

template<typename K, typename V>
void testThatAllElementsAreContained(const PerfectHashTable<K, V>& pht, const std::vector<std::pair<K,V>>& elements) {
  for ( const auto& element : elements ) {
    ASSERT_TRUE(pht.contains(element.first));
  }
}

template<typename K, typename V>
void testThatAllElementsAreNotContained(const PerfectHashTable<K, V>& pht, const std::vector<std::pair<K,V>>& elements) {
  for ( const auto& element : elements ) {
    ASSERT_FALSE(pht.contains(element.first));
  }
}

template<typename K, typename V>
void testThatAllValuesAreCorrect(const PerfectHashTable<K, V>& pht, const std::vector<std::pair<K,V>>& elements) {
  for ( const auto& element : elements ) {
    ASSERT_TRUE(pht.contains(element.first));
    ASSERT_EQ(*pht.get(element.first), element.second);
  }
}

TEST(APerfectHashTable, verifiesThatAllInsertedElementsAreContained) {
  std::vector<std::pair<int,int>> elements =
    { P(0, 1), P(4, 7), P(8, 1), P(11, 2), P(15, 21), P(42, 5), P(211, 23) };
  auto pht = construct(elements, 2048);
  testThatAllElementsAreContained(pht, elements);
}

TEST(APerfectHashTable, verifiesThatValuesAreCorrectlyStored) {
  std::vector<std::pair<int,int>> elements =
    { P(0, 1), P(4, 7), P(8, 1), P(11, 2), P(15, 21), P(42, 5), P(211, 23) };
  auto pht = construct(elements, 2048);
  testThatAllValuesAreCorrect(pht, elements);
}

TEST(APerfectHashTable, testSomeElementIfTheyAreNotContained) {
  std::vector<std::pair<int,int>> elements =
    { P(0, 1), P(4, 7), P(8, 1), P(11, 2), P(15, 21), P(42, 5), P(211, 23) };
  auto pht = construct(elements, 2048);

  std::vector<std::pair<int,int>> not_contained =
    { P(512, 0), P(612, 0), P(1203, 1), P(3, 0), P(5, 0), P(6, 0) };
  testThatAllElementsAreNotContained(pht, not_contained);
}

TEST(APerfectHashTable, testElementAfterReconstruct) {
  std::vector<std::pair<int,int>> elements_1 =
    { P(0, 1), P(4, 7), P(8, 1), P(11, 2), P(15, 21), P(42, 5), P(211, 23) };
  auto pht = construct(elements_1, 2048);

  std::vector<std::pair<int,int>> elements_2 =
    { P(512, 0), P(612, 0), P(1203, 1), P(3, 0), P(5, 0), P(6, 0) };
  construct(pht, elements_2);

  testThatAllElementsAreNotContained(pht, elements_1);
  testThatAllElementsAreContained(pht, elements_2);
}

TEST(APerfectHashTable, StressTest) {
  const int N = 10000;
  const float p = 0.5;
  std::vector<std::pair<int,int>> elements;
  for ( int i = 0; i < N; ++i ) {
    const float prob = utils::Randomize::instance().getRandomFloat(0.0f, 1.0f, sched_getcpu());
    if ( prob < p ) {
      elements.push_back(P(i, i));
    }
  }
  auto pht = construct(elements, N);
  testThatAllElementsAreContained(pht, elements);
  testThatAllValuesAreCorrect(pht, elements);
}


}  // namespace ds
}  // namespace mt_kahypar
