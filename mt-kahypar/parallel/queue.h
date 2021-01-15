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
class queue {
  public:
    queue(size_t max_size);
    queue(const queue& q);

    bool write(T data);

    bool read(T &data);

    void clear();

    bool deactivate();

  private:
    std::atomic<bool> writer_lock; // manages access to the writer queue
    std::atomic<bool> deactivated; // if the queue is deactivated, do not write to it
    std::vector<T> writer_queue;
    std::vector<T> reader_queue;
    size_t max_size; // max size of reader and writer queue to prevent reallocation
};

template<class T>
queue<T>::queue(size_t max_size) : max_size(max_size) {
  writer_queue.reserve(max_size);
  reader_queue.reserve(max_size);
  writer_lock = false;
}

template<class T>
queue<T>::queue(const queue& q) {
  writer_queue = q.writer_queue;
  reader_queue = q.reader_queue;
  writer_lock = false;
}

template<class T>
bool queue<T>::write(T data) {

  bool w_top = writer_lock.exchange(true, std::memory_order_acq_rel);
  bool deact = deactivated.load(std::memory_order_acq_rel);

  // if writer is locked or queue is deactivated, do not write to queue
  /* TODO: maybe separate to allow write retry if writer is locked <15-01-21, @noahares> */
  if (w_top || deact) {
    writer_lock.store(w_top, std::memory_order_release);
    return false;
  }

  // if writer queue is not full, append message to it
  if (writer_queue.size() < max_size) {
    writer_queue.push_back(data);
    writer_lock.store(w_top, std::memory_order_release);
    return true;

    /* TODO: ideally writer should be able to give reader its queue if it is full to start a new one,
       but then we would need a reader lock because the check for empty can occur while moving the writer to the reader queue.
       In this case acquiring locks in the same order should be kept in mind <15-01-21, @noahares> */
    /*
  } else if (reader_queue.empty()) {
    reader_queue = std::move(writer_queue);

    writer_queue = std::vector<T>();
    writer_queue.reserve(max_size);
    writer_queue.push_back(data); // append to queue
    writer_lock.store(w_top, std::memory_order_release);
    return true;
*/
  }
  writer_lock.store(w_top, std::memory_order_release); // restore writer's top
  return false;
}

template<class T>
bool queue<T>::read(T &data) {

  // if reader has nothing to read, try to acquire writer queue
  if (reader_queue.empty()) {
    bool w_top = writer_lock.exchange(true, std::memory_order_acq_rel);
    if (w_top || writer_queue.empty()) {
      writer_lock.store(w_top, std::memory_order_release);
      return false;
    } else {
      reader_queue = std::move(writer_queue);

      writer_queue = std::vector<T>();
      writer_queue.reserve(max_size);
      writer_lock.store(w_top, std::memory_order_release);
    }

  }
  data = reader_queue.back();
  reader_queue.pop_back();
  return true;
}

template<class T>
void queue<T>::clear() {
  // clean queues
  reader_queue.clear();
  writer_queue.clear();
  writer_lock.store(false, std::memory_order_acq_rel);
}

template<class T>
bool queue<T>::deactivate() {
  // deactivate the queue to prevent writes to it
  bool w_top = writer_lock.exchange(true, std::memory_order_acq_rel);
  if (w_top) {
    return false;
  }
  deactivated.store(true, std::memory_order_acq_rel);
  clear();
  return true;
}


} //namespace parallel
} //namespace kahypar
