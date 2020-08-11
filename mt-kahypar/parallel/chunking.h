#pragma once
#include <vector>

namespace chunking {
  template <typename T1, typename T2>
  auto idiv_ceil(T1 a, T2 b) {
    return static_cast<T1>((static_cast<unsigned long long>(a)+b-1) / b);
  };

  std::pair<size_t, size_t> bounds(size_t i, size_t n, size_t chunk_size) {
    return std::make_pair(i * chunk_size, std::min(n, (i+1) * chunk_size));
  };
}