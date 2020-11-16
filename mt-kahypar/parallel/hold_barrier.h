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

#include "mt-kahypar/macros.h"
#include <condition_variable>
#include <mutex>

namespace mt_kahypar {
namespace parallel {

/* HoldBarrier implements a two part barrier that requires _size threads to
 * enter it via aquire() and finally exit it via release(). _size can be reduced
 * via lowerSize() until all threads have entered. Having more than _size
 * threads call either aquire() or release() is not allowed. To reset the
 * barrier with a given size call reset(size), this should only be done if no
 * thread currently has entered the barrier but not jet exited it. */

class HoldBarrier {
public:
  HoldBarrier(size_t size) : _size(size) {}

  void aquire() {
    std::unique_lock<std::mutex> lock(_mutex);
    ASSERT(_current < _size);
    if (++_current == _size) {
      lock.unlock();
      _release = true;
      _cv.notify_all();
    } else {
      _cv.wait(lock, [this] { return _current == _size || _release; });
      lock.unlock();
    }
  }

  void release() {
    std::unique_lock<std::mutex> lock(_mutex);
    ASSERT(_current > 0);
    if (--_current == 0) {
      lock.unlock();
      _cv.notify_all();
    } else {
      _cv.wait(lock, [this] { return _current == 0; });
      lock.unlock();
    }
  }

  void lowerSize() {
    std::unique_lock<std::mutex> lock(_mutex);
    ASSERT(_size > 0);
    --_size;
    lock.unlock();
    _cv.notify_all();
  }
  void reset(size_t size) {
    ASSERT(_current == 0);
    _size = size;
  }

private:
  size_t _size;
  std::condition_variable _cv;
  std::mutex _mutex;
  size_t _current = 0;
  bool _release = false;
};

} // namespace parallel
} // namespace mt_kahypar
