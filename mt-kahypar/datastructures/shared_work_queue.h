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
#include "mt-kahypar/definitions.h"
#include "mt-kahypar/macros.h"
#include "mt-kahypar/parallel/atomic_wrapper.h"
#include "mt-kahypar/parallel/stl/scalable_vector.h"
#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/utils/range.h"
#include <algorithm>
#include <random>
#include <tbb/parallel_sort.h>
namespace mt_kahypar {

template <class T> struct SharedWorkQueue {
  vec<T> work_queue;
  CAtomic<size_t> front;
  CAtomic<size_t> back;
  size_t size;

  SharedWorkQueue(size_t size) :size(size) {
    work_queue.resize(size);
    front.store(0);
    back.store(0);
  }

  void clear() {
    work_queue.assign(size, 0);
    front.store(0);
    back.store(0);
  }

  void append(vec<T> &elements) {
    size_t pos = back.fetch_add(elements.size(), std::memory_order_acq_rel);
    for (size_t i = 0; i < elements.size(); ++i) {
      work_queue[pos + i] = elements[i];
    }
  }

  bool try_pop(T &dest) {
    size_t pos = front.fetch_add(1, std::memory_order_acq_rel);
    if (pos < back) {
      dest = work_queue[pos];
      return true;
    }
    return false;
  }

  std::optional<IteratorRange<typename vec<T>::iterator>> try_pop(size_t num_seeds) {
    size_t pos = front.fetch_add(num_seeds, std::memory_order_acq_rel);
    if (pos < back) {
      const auto begin = work_queue.begin() + pos;
      const auto end = std::min(begin + num_seeds, work_queue.begin() + back);
      return IteratorRange(begin, end);
    }
    return {};
  }

  size_t unsafe_size() const {
    return back.load(std::memory_order_relaxed) - front.load(std::memory_order_relaxed);
  }

  struct Comparator {
    const vec<Gain> &value_vector;

    Comparator(const vec<Gain> &val_vec) : value_vector(val_vec) {}

    bool operator()(const T& i1, const T& i2) {
      return value_vector[i1] > value_vector[i2];
    }
  };

  void shuffle(size_t seed) {
/*
    vec<int> randoms;
    vec<int> indices;
    randoms.resize(work_queue.size());
    indices.resize(work_queue.size());
*/
    std::mt19937 g(seed);
/*
    std::uniform_int_distribution<> distrib(0, std::numeric_limits<int>::max());
    auto get_rand = [&]() { return distrib(g); };
    std::generate(randoms.begin(), randoms.end(), get_rand);
    std::iota(indices.begin(), indices.end(), 0);
    tbb::parallel_sort(indices.begin(), indices.end(),
                       Comparator(randoms));
*/
    std::shuffle(work_queue.begin() + front, work_queue.begin() + back, g);
/*
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
*/
  }

  // very hacky thing
  void sortByGain(PartitionedHypergraph& phg, const Context& context) {
    vec<Gain> gains(phg.initialNumNodes());
    auto begin = work_queue.begin() + front;
    auto end = work_queue.begin() + back;
    for (auto i = begin; i < end; ++i) {
      HypernodeID u = *i;
      const HypernodeWeight wu = phg.nodeWeight(u);
      const PartitionID from = phg.partID(u);
      const HypernodeWeight from_weight = phg.partWeight(from);
      PartitionID to = kInvalidPartition;
      HyperedgeWeight to_penalty = std::numeric_limits<HyperedgeWeight>::max();
      HypernodeWeight best_to_weight = from_weight - wu;
      for (PartitionID i = 0; i < phg.k(); ++i) {
        if (i != from) {
          const HypernodeWeight to_weight = phg.partWeight(i);
          const HyperedgeWeight penalty = phg.moveToPenalty(u, i);
          if ( ( penalty < to_penalty || ( penalty == to_penalty && to_weight < best_to_weight ) ) &&
              to_weight + wu <= context.partition.max_part_weights[i] ) {
            to_penalty = penalty;
            to = i;
            best_to_weight = to_weight;
          }
        }
      }
      const Gain gain = to != kInvalidPartition ? phg.moveFromBenefit(u) - to_penalty : std::numeric_limits<HyperedgeWeight>::min();
      gains[u] = gain;
    }
    std::sort(work_queue.begin() + front, work_queue.begin() + back, Comparator(gains));
  }
};


} // namespace mt_kahypar
