/*******************************************************************************
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2019 Lars Gottesbüren <lars.gottesbueren@kit.edu>
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
#pragma once

#include <cstdint>
#include <functional>
#include <iostream>
#include <numeric>
#include <vector>
#include <cassert>

#include "mt-kahypar/parallel/chunking.h"
#include "mt-kahypar/parallel/parallel_prefix_sum.h"

namespace mt_kahypar::parallel {

// KeyFunc must be thread safe
// returns the bucket bounds
template <class InputRange, class OutputRange, class KeyFunc>
vec<uint32_t> counting_sort(const InputRange& input, OutputRange& output,
                            size_t max_num_buckets, KeyFunc& get_bucket, size_t num_tasks) {

  vec<uint32_t> global_bucket_begins(max_num_buckets + 2, 0);

  const size_t n = input.size();

  if (num_tasks > 1 && n > (1 << 17)) {
    const size_t chunk_size = chunking::idiv_ceil(n, num_tasks);

    // thread local counting
    vec<vec<uint32_t>> thread_local_bucket_ends(num_tasks);   // use vector of vector to avoid false sharing. maybe even task-local vector and then copy?
    tbb::parallel_for(size_t(0), num_tasks, [&](const size_t taskID) {
      vec<uint32_t>& bucket_ends = thread_local_bucket_ends[taskID];
      bucket_ends.resize(max_num_buckets, 0);
      for (auto[i,last] = chunking::bounds(taskID, n, chunk_size); i < last; ++i) {
        bucket_ends[get_bucket(input[i])]++;
      }
    });

    // prefix sum local bucket sizes for local offsets
    if (max_num_buckets > 1 << 10) {
      tbb::parallel_for(0UL, max_num_buckets, [&](size_t bucket) {
        for (size_t i = 1; i < num_tasks; ++i) {
          thread_local_bucket_ends[i][bucket] += thread_local_bucket_ends[i - 1][bucket]; // EVIL for locality!
        }
      });
    } else {
      for (size_t bucket = 0; bucket < max_num_buckets; ++bucket) {
        for (size_t i = 1; i < num_tasks; ++i) {
          thread_local_bucket_ends[i][bucket] += thread_local_bucket_ends[i - 1][bucket]; // EVIL for locality!
        }
      }
    }

    // prefix sum over bucket
    assert(global_bucket_begins.size()  >= thread_local_bucket_ends.back().size() + 1);
    parallel_prefix_sum(thread_local_bucket_ends.back().cbegin(), thread_local_bucket_ends.back().cend(),
                        global_bucket_begins.begin() + 1,
                        std::plus<>(), 0);

    // element assignment
    tbb::parallel_for(size_t(0), num_tasks, [&](const size_t taskID) {
      vec<uint32_t>& bucketEnds = thread_local_bucket_ends[taskID];
      // reverse iteration makes the algorithm stable
      for (auto [first,i] = chunking::bounds(taskID, n, chunk_size); i > first; --i) {
        size_t bucket = get_bucket(input[i-1]);
        output[global_bucket_begins[bucket] + (--bucketEnds[bucket])] = input[i-1];
      }
    });

  } else {
    for (size_t i = 0; i < input.size(); ++i) global_bucket_begins[get_bucket(input[i]) + 2]++;
    std::partial_sum(global_bucket_begins.begin(), global_bucket_begins.end(), global_bucket_begins.begin());
    for (size_t i = 0; i < input.size(); ++i) output[global_bucket_begins[get_bucket(input[i]) + 1]++] = input[i] ;
  }

  global_bucket_begins.pop_back();    // did the +2 trick
  return global_bucket_begins;
}

}  // namespace mt_kahypar::parallel
