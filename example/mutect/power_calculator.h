//
// Created by zp on 11/23/16.
//

#ifndef MUTECT_POWER_CALCULATOR_H_
#define MUTECT_POWER_CALCULATOR_H_

#include <cmath>
#include <functional>

namespace ncic {
namespace mutect {

struct PowerCacheKey {
  int n;
  double delta;

  PowerCacheKey(int n_, double delta_) {
    n = n_;
    delta = delta_;
  }

  bool operator==(const PowerCacheKey& rhs) const {
    if (&rhs == this) return true;
    const int EPS = 1e-8;
    if (n != rhs.n) return false;
    if (std::abs(delta - rhs.delta) > EPS) return false;
    return true;
  }
};

} // mutect
} // ncic

// the hash function for PowerCacheKey
namespace std {
template<> struct hash<ncic::mutect::PowerCacheKey>
{
  size_t operator()(ncic::mutect::PowerCacheKey const& rhs) const
  {
    size_t const h1(std::hash<int>{}(rhs.n));
    size_t const h2(std::hash<double>{}(rhs.delta));
    return h1 ^ h2;
  }
};
}

#endif
