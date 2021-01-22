/*******************************************************************************
 * This file is part of MT-KaHyPar.
 *
 * Copyright (C) 2021 Noah Wahl <noah.wahl@student.kit.edu>
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

#include "mt-kahypar/datastructures/hypergraph_common.h"
#include "mt-kahypar/macros.h"
#include "mt-kahypar/parallel/atomic_wrapper.h"
#include "mt-kahypar/parallel/stl/scalable_vector.h"
#include <algorithm>
#include <tbb/parallel_sort.h>
namespace mt_kahypar {

struct Comparator {
  const vec<int> &value_vector;

  Comparator(const vec<int> &val_vec) : value_vector(val_vec) {}

  bool operator()(int i1, int i2) {
    return value_vector[i1] < value_vector[i2];
  }
};

template <class T> struct SharedWorkQueue {
  vec<T> work_queue;
  CAtomic<size_t> front;

  SharedWorkQueue() { front.store(0); }

  void clear() {
    work_queue.clear();
    front.store(0);
  }

  void append(vec<T> &elements) {
    work_queue.insert(work_queue.end(), elements.begin(), elements.end());
  }

  bool try_pop(T &dest) {
    size_t pos = front.fetch_add(1, std::memory_order_acq_rel);
    if (pos < work_queue.size()) {
      dest = work_queue[pos];
      return true;
    }
    return false;
  }

  bool try_pop(vec<T> &dest, size_t num_seeds) {
    size_t pos = front.fetch_add(num_seeds, std::memory_order_acq_rel);
    if (pos < work_queue.size()) {
      for (size_t i = 0; (i < num_seeds) && (pos + i < work_queue.size());
           ++i) {
        dest.push_back(work_queue[pos + i]);
      }
      return true;
    }
    return false;
  }

  size_t unsafe_size() const {
    return work_queue.size() - front.load(std::memory_order_relaxed);
  }

  void shuffle() {
    vec<int> randoms;
    randoms.reserve(work_queue.size());
    std::srand(time(0));
    std::generate(randoms.begin(), randoms.end(), std::rand);
    tbb::parallel_sort(work_queue.begin(), work_queue.end(),
                       Comparator(randoms));
    tbb::parallel_for(tbb::blocked_range<HypernodeID>(0, work_queue.size()),
                      [&](const tbb::blocked_range<HypernodeID> &r) {
                        for (size_t i = r.begin(); i < r.end(); ++i) {
                          size_t swap =
                              r.begin() +
                              (std::rand() %
                               static_cast<size_t>(r.end() - r.begin()));
                          ASSERT(swap >= r.begin() && swap < r.end());
                          std::swap(work_queue[i], work_queue[swap]);
                        }
                      });
  }
};

} // namespace mt_kahypar
