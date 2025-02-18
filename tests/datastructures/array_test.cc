/*******************************************************************************
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2019 Tobias Heuer <tobias.heuer@kit.edu>
 *
 * Mt-KaHyPar is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Mt-KaHyPar is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Mt-KaHyPar.  If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/

#include <numeric>
#include <algorithm>

#include "gmock/gmock.h"

#include "mt-kahypar/parallel/stl/scalable_vector.h"
#include "mt-kahypar/datastructures/array.h"
#include "mt-kahypar/parallel/hardware_topology.h"
#include "mt-kahypar/parallel/tbb_initializer.h"
#include "tests/parallel/topology_mock.h"

using ::testing::Test;

namespace mt_kahypar {
namespace ds {

using TopoMock = mt_kahypar::parallel::TopologyMock<2>;
using HwTopology = mt_kahypar::parallel::HardwareTopology<TopoMock, parallel::topology_t, parallel::node_t>;
using TBB = mt_kahypar::parallel::TBBInitializer<HwTopology, false>;

TEST(AArray, WritesAnValueToStrippedVector1) {
  Array<int> vec(256, 0);
  vec[0] = 31;
  ASSERT_EQ(31, vec[0]);
}

TEST(AArray, WritesAnValueToStrippedVector2) {
  Array<int> vec(256, 0);
  vec[65] = 35;
  ASSERT_EQ(35, vec[65]);
}

TEST(AArray, WritesAnValueToStrippedVector3) {
  Array<int> vec(256, 0);
  vec[127] = 42;
  ASSERT_EQ(42, vec[127]);
}

TEST(AArray, WritesAnValueToStrippedVector4) {
  Array<int> vec(256, 0);
  vec[128] = 43;
  ASSERT_EQ(43, vec[128]);
}

TEST(AArray, WritesValuesToWholeVector) {
  Array<int> vec(256, 0);

  for ( size_t i = 0; i < vec.size(); ++i ) {
    vec[i] = i;
  }
  for ( size_t i = 0; i < vec.size(); ++i ) {
    ASSERT_EQ(i, vec[i]);
  }
}

TEST(AArray, IsInitializedWithNonDefaultValues) {
  Array<int> vec(256, 42);

  for ( size_t i = 0; i < vec.size(); ++i ) {
    ASSERT_EQ(42, vec[i]);
  }
}

TEST(AArray, AssignValuesToAlreadyInitializedVector1) {
  Array<int> vec(256, 0);

  const size_t count = 31;
  vec.assign(count, 42);
  for ( size_t i = 0; i < count; ++i ) {
    ASSERT_EQ(42, vec[i]);
  }
  for ( size_t i = count; i < vec.size(); ++i ) {
    ASSERT_EQ(0, vec[i]);
  }
}

TEST(AArray, AssignValuesToAlreadyInitializedVector2) {
  Array<int> vec(256, 0);

  const size_t count = 42;
  vec.assign(count, 42);
  for ( size_t i = 0; i < count; ++i ) {
    ASSERT_EQ(42, vec[i]);
  }
  for ( size_t i = count; i < vec.size(); ++i ) {
    ASSERT_EQ(0, vec[i]);
  }
}

TEST(AArray, AssignValuesToAlreadyInitializedVector3) {
  Array<int> vec(256, 0);

  const size_t count = 127;
  vec.assign(count, 42);
  for ( size_t i = 0; i < count; ++i ) {
    ASSERT_EQ(42, vec[i]);
  }
  for ( size_t i = count; i < vec.size(); ++i ) {
    ASSERT_EQ(0, vec[i]);
  }
}

TEST(AArray, AssignValuesToAlreadyInitializedVector4) {
  Array<int> vec(256, 0);

  const size_t count = 128;
  vec.assign(count, 42);
  for ( size_t i = 0; i < count; ++i ) {
    ASSERT_EQ(42, vec[i]);
  }
  for ( size_t i = count; i < vec.size(); ++i ) {
    ASSERT_EQ(0, vec[i]);
  }
}

TEST(AArray, AssignValuesToAlreadyInitializedVector5) {
  Array<int> vec(256, 0);

  const size_t count = 256;
  vec.assign(count, 42);
  for ( size_t i = 0; i < count; ++i ) {
    ASSERT_EQ(42, vec[i]);
  }
  for ( size_t i = count; i < vec.size(); ++i ) {
    ASSERT_EQ(0, vec[i]);
  }
}

TEST(AArray, FilledWithNumbersFromZeroToN) {
  Array<int> vec(256, 0);
  std::iota(vec.begin(), vec.end(), 0);
  for ( size_t i = 0; i < vec.size(); ++i ) {
    ASSERT_EQ(i, vec[i]);
  }
}

TEST(AArray, ChecksDistanceBetweenTwoPointers1) {
  Array<int> vec(256, 0);
  ASSERT_EQ(5, std::distance(vec.begin(), vec.begin() + 5));
}

TEST(AArray, ChecksDistanceBetweenTwoPointers2) {
  Array<int> vec(256, 0);
  ASSERT_EQ(42, std::distance(vec.begin() + 24, vec.begin() + 66));
}

TEST(AArray, ChecksDistanceBetweenTwoPointers3) {
  Array<int> vec(256, 0);
  ASSERT_EQ(256, std::distance(vec.begin(), vec.end()));
}

TEST(AArray, CanBeSorted) {
  Array<int> vec(256, 0);
  for ( size_t i = 0; i < vec.size(); ++i ) {
    vec[i] = (vec.size() - 1) - i;
  }
  std::sort(vec.begin(), vec.end());
  for ( size_t i = 0; i < vec.size(); ++i ) {
    ASSERT_EQ(i, vec[i]);
  }
}



TEST(AArray, MemcopiesContentToVector) {
  Array<int> vec(256, 0);
  std::vector<int> vec2(256, 5);
  memcpy(vec.data(), vec2.data(), sizeof(int) * 256);
  for ( size_t i = 0; i < vec.size(); ++i ) {
    ASSERT_EQ(5, vec[i]);
  }
}

TEST(AArray, IsInitializedWithMemoryChunkFromMemoryPool) {
  parallel::MemoryPool::instance().deactivate_minimum_allocation_size();
  parallel::MemoryPool::instance().register_memory_group("TEST_GROUP", 1);
  parallel::MemoryPool::instance().register_memory_chunk("TEST_GROUP", "TEST_CHUNK", 5, sizeof(size_t));
  parallel::MemoryPool::instance().allocate_memory_chunks();

  Array<size_t> vec("TEST_GROUP", "TEST_CHUNK", 5);
  ASSERT_EQ(parallel::MemoryPool::instance().mem_chunk("TEST_GROUP", "TEST_CHUNK"), (char *) vec.data());

  parallel::MemoryPool::instance().free_memory_chunks();
}

TEST(AArray, IsInitializedWithSeveralMemoryChunksFromMemoryPool) {
  parallel::MemoryPool::instance().deactivate_minimum_allocation_size();
  parallel::MemoryPool::instance().register_memory_group("TEST_GROUP", 1);
  parallel::MemoryPool::instance().register_memory_chunk("TEST_GROUP", "TEST_CHUNK_1", 5, sizeof(size_t));
  parallel::MemoryPool::instance().register_memory_chunk("TEST_GROUP", "TEST_CHUNK_2", 5, sizeof(size_t));
  parallel::MemoryPool::instance().allocate_memory_chunks();

  Array<size_t> vec_1("TEST_GROUP", "TEST_CHUNK_1", 5);
  Array<size_t> vec_2("TEST_GROUP", "TEST_CHUNK_2", 5);
  ASSERT_EQ(parallel::MemoryPool::instance().mem_chunk("TEST_GROUP", "TEST_CHUNK_1"), (char *) vec_1.data());
  ASSERT_EQ(parallel::MemoryPool::instance().mem_chunk("TEST_GROUP", "TEST_CHUNK_2"), (char *) vec_2.data());
  ASSERT_NE(vec_1.data(), vec_2.data());

  parallel::MemoryPool::instance().free_memory_chunks();
}

TEST(AArray, ReleasesMemoryChunkFromMemoryPoolInDestructor) {
  parallel::MemoryPool::instance().deactivate_minimum_allocation_size();
  parallel::MemoryPool::instance().register_memory_group("TEST_GROUP", 1);
  parallel::MemoryPool::instance().register_memory_chunk("TEST_GROUP", "TEST_CHUNK", 5, sizeof(size_t));
  parallel::MemoryPool::instance().allocate_memory_chunks();

  {
    Array<size_t> vec("TEST_GROUP", "TEST_CHUNK", 5);
    ASSERT_EQ(parallel::MemoryPool::instance().mem_chunk("TEST_GROUP", "TEST_CHUNK"), (char *) vec.data());
    ASSERT_EQ(nullptr, parallel::MemoryPool::instance().request_mem_chunk("TEST_GROUP", "TEST_CHUNK", 5, sizeof(size_t)));
  }

  ASSERT_NE(nullptr, parallel::MemoryPool::instance().request_mem_chunk("TEST_GROUP", "TEST_CHUNK", 5, sizeof(size_t)));

  parallel::MemoryPool::instance().free_memory_chunks();
}

TEST(AArray, AllocatesOwnMemoryIfNotAvailableInMemoryPool) {
  parallel::MemoryPool::instance().deactivate_minimum_allocation_size();
  parallel::MemoryPool::instance().register_memory_group("TEST_GROUP", 1);
  parallel::MemoryPool::instance().register_memory_chunk("TEST_GROUP", "TEST_CHUNK", 5, sizeof(size_t));
  parallel::MemoryPool::instance().allocate_memory_chunks();

  Array<size_t> vec("TEST_GROUP", "OTHER_CHUNK", 5);
  ASSERT_NE(nullptr, vec.data());
  ASSERT_NE(parallel::MemoryPool::instance().mem_chunk("TEST_GROUP", "TEST_CHUNK"), (char *) vec.data());

  parallel::MemoryPool::instance().free_memory_chunks();
}

TEST(AArray, AllocatesOwnMemoryOnOverAllocationInMemoryPool) {
  parallel::MemoryPool::instance().deactivate_minimum_allocation_size();
  parallel::MemoryPool::instance().register_memory_group("TEST_GROUP", 1);
  parallel::MemoryPool::instance().register_memory_chunk("TEST_GROUP", "TEST_CHUNK", 5, sizeof(size_t));
  parallel::MemoryPool::instance().allocate_memory_chunks();

  Array<size_t> vec("TEST_GROUP", "TEST_CHUNK", 10);
  ASSERT_NE(nullptr, vec.data());
  ASSERT_NE(parallel::MemoryPool::instance().mem_chunk("TEST_GROUP", "TEST_CHUNK"), (char *) vec.data());

  parallel::MemoryPool::instance().free_memory_chunks();
}

TEST(AArray, AllocatesOwnMemoryIfAlreadyRequestedInMemoryPool) {
  parallel::MemoryPool::instance().deactivate_minimum_allocation_size();
  parallel::MemoryPool::instance().register_memory_group("TEST_GROUP", 1);
  parallel::MemoryPool::instance().register_memory_chunk("TEST_GROUP", "TEST_CHUNK", 5, sizeof(size_t));
  parallel::MemoryPool::instance().allocate_memory_chunks();

  Array<size_t> vec_1("TEST_GROUP", "TEST_CHUNK", 5);
  Array<size_t> vec_2("TEST_GROUP", "TEST_CHUNK", 5);
  ASSERT_EQ(parallel::MemoryPool::instance().mem_chunk("TEST_GROUP", "TEST_CHUNK"), (char *) vec_1.data());
  ASSERT_NE(nullptr, vec_2.data());
  ASSERT_NE(parallel::MemoryPool::instance().mem_chunk("TEST_GROUP", "TEST_CHUNK"), (char *) vec_2.data());

  parallel::MemoryPool::instance().free_memory_chunks();
}



}  // namespace ds
}  // namespace mt_kahypar
