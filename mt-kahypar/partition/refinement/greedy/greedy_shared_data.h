/*******************************************************************************
 * This file is part of MT-KaHyPar.
 *
 * Copyright (C) 2020 Lars Gottesbüren <lars.gottesbueren@kit.edu>
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

#include "mt-kahypar/datastructures/hypergraph_common.h"
#include <condition_variable>
#include <mutex>
namespace mt_kahypar {

class HoldBarrier {
public:
  HoldBarrier(size_t size) : _size(size) {}
  bool aquire() {
    std::unique_lock<std::mutex> lock{_mutex};
    if (_current == _size - 1) {
      ++_current;
      _cv.notify_all();
    } else if (_current == _size) {
      lock.unlock();
      return false;
    } else {
      ++_current;
      _cv.wait(lock, [this] { return _current == _size; });
    }
    return true;
  }
  bool release() {
    std::unique_lock<std::mutex> lock{_mutex};
    if (_current == 1) {
      --_current;
      _cv.notify_all();
    } else if (_current == 0) {
      lock.unlock();
      return false;
    } else {
      --_current;
      _cv.wait(lock, [this] { return _current == 0; });
    }
    return true;
  }
  bool lowerSize() {
    std::unique_lock<std::mutex> lock{_mutex};
    if (_size > 1) {
      --_size;
      _cv.notify_all();
      lock.unlock();
      return true;
    } else {
      lock.unlock();
      return false;
    }
  }
  void reset(size_t size) { _size = size; }

private:
  size_t _size;
  std::condition_variable _cv;
  std::mutex _mutex;
  size_t _current = 0;
};
using HypernodeIDMessageMatrix = vec<vec<HypernodeID>>;
struct GreedySharedData {
  HypernodeIDMessageMatrix messages;
  HoldBarrier hold_barrier;
  GreedySharedData(size_t num_threads)
      : messages(num_threads * num_threads, vec<HypernodeID>()),
        hold_barrier(num_threads) {}
};
} // namespace mt_kahypar
