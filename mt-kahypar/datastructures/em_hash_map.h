/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2019 Tobias Heuer <tobias.heuer@kit.edu>
 * Copyright (C) 2019 Lars Gottesbüren <lars.gottesbueren@kit.edu>
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

#include "ehash_table.h"

#include "kahypar/meta/mandatory.h"

#include "mt-kahypar/macros.h"

namespace mt_kahypar {
namespace ds {

template <typename Key = Mandatory>
struct chash {
    MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE static size_t hash_f(size_t x) {
        return x ^ (x << 7);
        // x += 0x9e3779b97f4a7c15;
        // x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9;
        // // x = (x ^ (x >> 27)) * 0x94d049bb133111eb;
        // return x ^ (x >> 31);
    }
    MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE size_t operator()(Key x) const { return hash_f(x); }
};

template <typename Key = Mandatory,
          typename Value = Mandatory,
          typename Hash = std::hash<Key>>
class EmHashMap {
    using HashMap = emhash7::HashMap<Key, Value, Hash>;

 public:
  explicit EmHashMap() :
    _map(4, 0.85) {
      // TODO: allocate large chunk?
  }

  EmHashMap(const EmHashMap&) = delete;
  EmHashMap& operator= (const EmHashMap& other) = delete;

  EmHashMap(EmHashMap&& other) :
    _map(std::move(other._map)) {
  }

  ~EmHashMap() = default;

//   size_t capacity() const {
//       // TODO
//   }

  size_t size() const {
    return _map.size();
  }

  typename HashMap::const_iterator begin() const {
    return _map.cbegin();
  }

  typename HashMap::const_iterator end() const {
    return _map.cend();
  }

  typename HashMap::iterator begin() {
    return _map.begin();
  }

  typename HashMap::iterator end() {
    return _map.end();
  }

  bool contains(const Key key) const {
    return _map.contains(key);
  }

  Value& operator[] (const Key key) {
    Value* value = _map.try_get(key);
    if (value == nullptr) {
      auto result = _map.insert(key, Value());
      return result.first->second;
    } else {
      return *value;
    }
  }

  const Value & get(const Key key) const {
    ASSERT(contains(key));
    return _map[key];
  }

  const Value* get_if_contained(const Key key) const {
    return _map.try_get(key);
  }

  void clear() {
    _map.clear();
  }

  void freeInternalData() {
    // TODO: ?
    _map.clear();
  }

  size_t size_in_bytes() const {
    return _map.size_in_bytes();
  }

 private:
  HashMap _map;
};

} // namespace ds
} // namespace mt_kahypar
