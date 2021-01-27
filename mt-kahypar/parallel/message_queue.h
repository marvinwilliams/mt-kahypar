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

#include <atomic>
#include <queue>

namespace mt_kahypar {
namespace parallel {

template<class T>
class message_queue {
  public:
    message_queue();
    message_queue(const message_queue& q);

    bool try_write(T data);

    bool read(T &data);

    void unsafe_clear();

    bool try_clear();

  private:
    std::atomic<bool> writer_lock; // manages access to the writer queue
    std::vector<T> writer_queue;
    std::vector<T> reader_queue;
};

template<class T>
message_queue<T>::message_queue() {
  writer_lock = false;
}

template<class T>
message_queue<T>::message_queue(const message_queue& q) {
  writer_queue = q.writer_queue;
  reader_queue = q.reader_queue;
  writer_lock = false;
}

template<class T>
bool message_queue<T>::try_write(T data) {

  bool w_top = writer_lock.exchange(true, std::memory_order_acq_rel);

  // if writer is locked or queue is deactivated, do not write to queue
  if (w_top) {
    return false;
  }

  writer_queue.push_back(data);
  writer_lock.store(w_top, std::memory_order_release);
  return true;
}

template<class T>
bool message_queue<T>::read(T &data) {

  // if reader has nothing to read, try to acquire writer queue
  if (reader_queue.empty()) {
    bool w_top = writer_lock.exchange(true, std::memory_order_acq_rel);
    if (w_top) {
      return false;
    }
    if (writer_queue.empty()) {
      writer_lock.store(w_top, std::memory_order_release);
      return false;
    } else {
      reader_queue.swap(writer_queue);
      writer_lock.store(w_top, std::memory_order_release);
    }

  }
  data = reader_queue.back();
  reader_queue.pop_back();
  return true;
}

template<class T>
void message_queue<T>::unsafe_clear() {
  // clean queues
  reader_queue.clear();
  writer_queue.clear();
}

template<class T>
bool message_queue<T>::try_clear() {
  // deactivate the queue to prevent writes to it
  bool w_top = writer_lock.exchange(true, std::memory_order_acquire);
  if (w_top) {
    return false;
  }
  unsafe_clear();
  writer_lock.store(w_top, std::memory_order_release);
  return true;
}

} //namespace parallel
} //namespace kahypar
