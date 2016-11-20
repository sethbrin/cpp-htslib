//
// Created by zp on 9/3/16.
//

#include "easehts/utils.h"

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

} // util
} // easehts
} // ncic
