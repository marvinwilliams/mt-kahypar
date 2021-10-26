/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2019 Lars Gottesbüren <lars.gottesbueren@kit.edu>
 * Copyright (C) 2020 Tobias Heuer <tobias.heuer@kit.edu>
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

#include <type_traits>

#include "tbb/concurrent_vector.h"

#include "mt-kahypar/parallel/stl/scalable_vector.h"

namespace mt_kahypar {
namespace ds {

template<typename Key, typename Value>
class PerfectHashTable {

  static constexpr bool debug = false;
  static constexpr bool enable_heavy_assert = false;

  static std::vector<size_t> primes;

  struct Element {
    Key key;
    Value value;
    uint16_t threshold = 0;
  };

  struct Bucket {
    explicit Bucket() :
      start(0),
      size(0),
      hash(0),
      threshold(0) { }

    size_t numRequiredSlots() const {
      return size * size;
    }

    size_t start;
    size_t size;
    Key hash;
    uint16_t threshold;
  };

  static_assert(std::is_integral<Key>(), "Key is not an integral type");
  static_assert(std::is_trivially_copyable<Element>(), "Elements are not trivially copyable");

 public:
  explicit PerfectHashTable(const Key universe_size) :
    _is_constructed(false),
    _tmp_elements(),
    _threshold(1),
    _prime(0),
    _first_level_hash(0),
    _buckets(),
    _slots(),
    _tmp_bucket_elements() {
    if (static_cast<size_t>(universe_size) >= primes.back()) {
      ERROR("Perfect Hash Table can not handle such large universe size");
    }

    for ( size_t i = 0; i < primes.size(); ++i ) {
      if ( primes[i] > static_cast<size_t>(universe_size) ) {
        _prime = primes[i];
        break;
      }
    }
  }

  void add(const Key key, const Value value) {
    ASSERT(!_is_constructed);
    _tmp_elements.emplace_back(Element { key, value, 0 });
  }

  const Value* get(const Key key) const {
    ASSERT(_is_constructed);
    const size_t idx = firstLevelHash(key);
    if ( !isEmpty(idx) ) {
      const Bucket& bucket = _buckets[idx];
      const size_t pos = bucket.start + secondLevelHash(key, bucket);
      return ( _slots[pos].threshold == _threshold &&
        _slots[pos].key == key ) ? &_slots[pos].value : nullptr;
    } else {
      return nullptr;
    }
  }

  bool contains(const Key key) const {
    ASSERT(_is_constructed);
    const size_t idx = firstLevelHash(key);
    if ( !isEmpty(idx) ) {
      const Bucket& bucket = _buckets[idx];
      const size_t pos = bucket.start + secondLevelHash(key, bucket);
      return _slots[pos].threshold == _threshold && _slots[pos].key == key;
    } else {
      return false;
    }
  }

  void clear() {
    _is_constructed = false;
    _tmp_elements.clear();
    clearBuckets();
  }

  void construct() {
    if ( _buckets.size() < _tmp_elements.size() ) {
      _buckets.resize(_tmp_elements.size());
    }

    DBG << "Construct first-level hash function";
    const size_t n = _tmp_elements.size();
    _first_level_hash = 3;
    size_t num_required_slots = testFirstLevelHash();
    while ( num_required_slots > 5 * n ) {
      ++_first_level_hash;
      num_required_slots = testFirstLevelHash();
    }
    DBG << "First-Level Hash =" << _first_level_hash
        << ", Num Required Slots =" << num_required_slots
        << ", Max Allowed Slots =" << (5*n);

    if ( _slots.size() < num_required_slots ) {
      _slots.resize(num_required_slots);
    }

    size_t start = 0;
    for ( size_t i = 0; i < n; ++i ) {
      if ( _buckets[i].threshold == _threshold ) {
        _buckets[i].start = start;
        start += _buckets[i].numRequiredSlots();
      }
    }

    DBG << "Construct second-level hash functions";
    for ( const Element& element : _tmp_elements ) {
      const size_t bucket = firstLevelHash(element.key);
      ASSERT(!isEmpty(bucket));
      _slots[_buckets[bucket].start++] = element;
    }

    for ( size_t i = 0; i < n; ++i ) {
      if ( !isEmpty(i) ) {
        _buckets[i].start -= _buckets[i].size;
        constructSecondLevelHashFunction(i);
      }
    }

    _is_constructed = true;
  }

  size_t size() const {
    return _tmp_elements.size();
  }

 private:
  size_t testFirstLevelHash() {
    DBG << "Test first-level hash =" << _first_level_hash;

    clearBuckets();
    const size_t n = _tmp_elements.size();
    size_t num_required_slots = 0;
    for ( const Element& element : _tmp_elements ) {
      const size_t i = static_cast<size_t>(_first_level_hash * element.key) % n;
      if ( _buckets[i].threshold < _threshold ) {
        _buckets[i].threshold = _threshold;
        _buckets[i].size = 0;
      }
      num_required_slots -= _buckets[i].numRequiredSlots();
      ++_buckets[i].size;
      num_required_slots += _buckets[i].numRequiredSlots();
    }

    return num_required_slots;
  }

  bool testSecondLevelHashFunction(const Bucket& bucket) {
    const size_t start = bucket.start;
    for ( size_t i = 0; i < _tmp_bucket_elements.size(); ++i ) {
      const Element& elem = _tmp_bucket_elements[i];
      const size_t pos = start + secondLevelHash(elem.key, bucket);
      ASSERT(pos < start + bucket.numRequiredSlots());
      if ( _slots[pos].threshold < _threshold  ) {
        _slots[pos] = Element { elem.key, elem.value, _threshold };
      } else /* collision */ {
        for ( int j = i - 1; j >= 0; --j ) {
          const Element& e = _tmp_bucket_elements[j];
          const size_t pos = start + secondLevelHash(e.key, bucket);
          _slots[pos].threshold = 0;
        }
        return false;
      }
    }
    return true;
  }

  void constructSecondLevelHashFunction(const size_t i) {
    ASSERT(!isEmpty(i));
    DBG << "Construct second-level hash function in bucket" << i;
    _tmp_bucket_elements.clear();
    for ( size_t j = 0; j < _buckets[i].size; ++j ) {
      _tmp_bucket_elements.push_back(_slots[_buckets[i].start + j]);
    }

    _buckets[i].hash = 3;
    while ( !testSecondLevelHashFunction(_buckets[i]) ) {
      ++_buckets[i].hash;
    }

    DBG << "Second-level hash for bucket" << i << "is" << _buckets[i].hash;
  }

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE size_t firstLevelHash(const Key key) const {
    return ( static_cast<size_t>(_first_level_hash) *
      static_cast<size_t>(key) ) % _tmp_elements.size();
  }

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE size_t secondLevelHash(const Key key,
                                                            const Bucket& bucket) const {
    return ( ( static_cast<size_t>(bucket.hash) *
      static_cast<size_t>(key) ) % _prime ) % bucket.numRequiredSlots();
  }

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE bool isEmpty(const size_t i) const {
    ASSERT(i < _buckets.size());
    return _buckets[i].threshold < _threshold;
  }

  void clearBuckets() {
    ++_threshold;
    if ( _threshold == std::numeric_limits<uint16_t>::max() ) {
      _threshold = 1;
      for ( size_t i = 0; i < _buckets.size(); ++i ) {
        _buckets[i].threshold = 0;
      }
    }
  }

  bool _is_constructed;
  tbb::concurrent_vector<Element> _tmp_elements;

  uint16_t _threshold;
  size_t _prime;
  Key _first_level_hash;
  vec<Bucket> _buckets;
  vec<Element> _slots;
  vec<Element> _tmp_bucket_elements;
};

template<typename Key, typename Value>
std::vector<size_t> PerfectHashTable<Key, Value>::primes =
  { 11, 19, 37, 47, 53, 59, 67, 73, 79, 83,
    157, 347, 449, 661, 701, 827, 829, 937, 971, 977,
    1223, 3607, 3881, 4013, 4519, 6553, 7487, 9521, 9749, 9901,
    17599, 20023, 20353, 23741, 31123, 38833, 41479, 42979, 59581, 88997,
    163393, 199931, 453707, 574699, 639053, 761711, 769151, 770897, 854897, 889907,
    1145983, 1551019, 2013289, 2248291, 2720009, 3747181, 4997771, 6153109, 7931141, 8012839,
    13519853, 17438917, 33453223, 36783401, 43548833, 50692867, 52247399, 52853393, 72598349, 99975103, 105701527 };

}  // namespace ds
}  // namespace mt_kahypar
