/**
 * Lock free queue for 1 writer and 1 reader threads.
 *
 * Writer and reader are working with a separate queues.
 * When writer writes to its own queue it checks if reader has anything to read.
 * If not then writer pass its queue to reader and starts new one for itself.
 *
 * The only place were reader and writer touch each other is when writer gives
 * its queue to reader via readerTop.
 */
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

  private:
    std::atomic<T*> writer_ptr;
    std::atomic<T*> reader_ptr;
    std::vector<T> writer_queue;
    std::vector<T> reader_queue;
    size_t max_size;
};

template<class T>
queue<T>::queue(size_t max_size) : max_size(max_size) {
  writer_queue.reserve(max_size);
  reader_queue.reserve(max_size);
  writer_ptr = &writer_queue.front();
  reader_ptr = &reader_queue.front();
}

template<class T>
queue<T>::queue(const queue& q) {
  writer_queue = q.writer_queue;
  reader_queue = q.reader_queue;
  writer_ptr = &writer_queue.front();
  reader_ptr = &reader_queue.front();
}

/*
 * Write data to the queue.
 * Algorithm:
 * 1. Retrieve writer top using atomic::exchange(null).
 *    This prevents reader from trying to take ownership of writers subqueue.
 * 2. If it is null then create new item.
 * 3. Otherwise add data to the end.
 * 4. Retrieve reader top using atomic::load(null).
 *    Using load instead of exchange prevents blocking of reader's subqueue.
 * 5. If it is null then set it to the writer top.
 * 6. Otherwise restore writer's top.
 */
template<class T>
bool queue<T>::write(T data) {

  T* w_top = writer_ptr.exchange(nullptr, std::memory_order_acq_rel);

  if (w_top == nullptr) {
    return false;

  }
  if (writer_queue.size() < max_size) {
    writer_queue.push_back(data); // append to queue
    writer_ptr.store(w_top, std::memory_order_release); // restore writer's top
    return true;
  } else if (reader_queue.empty()) { // reader don't have anything to read
    reader_queue = std::move(writer_queue);
    writer_queue = std::vector<T>();
    writer_queue.reserve(max_size);
    reader_ptr.store(w_top, std::memory_order_release); // give reader writer's queue

    writer_queue.push_back(data); // append to queue
    writer_ptr.store(w_top, std::memory_order_release); // restore writer's top
    return true;
  }
  return false;
}

/*
 * Read data from the queue.
 * Algorithm:
 * 1. Retrieve reader top using atomic::load().
 *    Using load instead of exchange prevets writer queue from overwriting readers one while reader is working with it.
 * 2. If it is null then:
 * 2.1. Retrieve writer top using atomic::exchange(null).
 *      Using exchange garantees that only writer or reader is owning writer's queue at each moment of time.
 * 2.2. If it is not null then assign it to the reader top otherwise exit.
 * 3. Shift reader's top to the next and return original top data.
 */
template<class T>
bool queue<T>::read(T &data)
{
  T* r_top = reader_ptr.load(std::memory_order_acquire);
  if (reader_queue.empty()) {
    T* w_top = writer_ptr.exchange(nullptr, std::memory_order_acq_rel);
    if (w_top == nullptr || writer_queue.empty()) {
      return false;
    } else {
      reader_queue = std::move(writer_queue);
      reader_ptr.store(w_top, std::memory_order_release);
      writer_queue = std::vector<T>();
      writer_queue.reserve(max_size);
      w_top = &writer_queue.front();
      writer_ptr.store(w_top, std::memory_order_release);
    }

  }
  reader_ptr.store(++r_top, std::memory_order_release);
  data = *r_top;
  return true;
}

template<class T>
void queue<T>::clear() {
  // clean queues
  reader_queue.clear();
  writer_queue.clear();
  reader_ptr = &reader_queue.front();
  writer_ptr = &writer_queue.front();
}


} //namespace parallel
} //namespace kahypar
