//
// Created by zp on 9/3/16.
//

#include "easehts/utils.h"

#include <algorithm>
#include <random>

namespace ncic {

namespace easehts {

namespace utils {

void tokenize(const std::string &s, char c, std::vector<std::string> *res) {
  auto end = s.end();
  auto start = end;

  for (auto it = s.begin(); it != end; ++it) {
    if (*it != c) {
      if (start == end)
        start = it;
      continue;
    }
    if (start != end) {
      res->emplace_back(start, it);
      start = end;
    }
  }
  if (start != end)
    res->emplace_back(start, end);
}

std::vector<int> SampleIndicesWithoutReplacemement(int n, int k,
                                                   ThreadLocalRandom& rnd) {
  std::vector<int> chosen_balls;
  chosen_balls.reserve(n);
  for (int i = 0; i < n; i++) {
    chosen_balls.push_back(i);
  }

  Shuffle(chosen_balls, rnd);
  return std::vector<int>(chosen_balls.begin(), chosen_balls.begin() + k);
}

// Random
const long long ThreadLocalRandom::kMultiplier = 0x5DEECE66DLL;
const long long ThreadLocalRandom::kAddend = 0xBLL;
const long long ThreadLocalRandom::kMask = (1LL << 48) - 1;


std::string RoundNearestFormat(double val, size_t max_digits) {
  double EPS = std::pow(10.0, -1.0 * (max_digits + 1)) * 5;

  // if val equals to -0.0000000000003, then just print 0
  // while the following print that -0
  if (std::fabs(val) < EPS) return "0";

  char str[80];
  for (int i = 0; i < max_digits; i++) {
    double round_val = RoundNearest(val, i);
    if (std::fabs(round_val - val) < EPS) {
      std::sprintf(str, ("%." + std::to_string(i) + "lf").c_str(), round_val);
      return std::string(str);
    }
  }
  std::sprintf(str, ("%." + std::to_string(max_digits) + "lf").c_str(), val);
  return std::string(str);
}

} // util
} // easehts
} // ncic
