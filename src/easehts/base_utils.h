//
// Created by zp on 11/25/16.
//

#ifndef EASEHTSLIB_BASE_UTILS_H_
#define EASEHTSLIB_BASE_UTILS_H_

#include "utils.h"

#include <vector>

namespace ncic {
namespace easehts {

/**
 * BaseUtils contains some basic utilities for manipulating nucleotides
 */
class BaseUtils {
 public:
  /**
   * Convert a simple base to a base index
   *
   * @param base [AaCcGgTt]
   * @return 0,1,2,3 or -1 if the base can't be understood
   */
  static int SimpleBaseToBaseIndex(char base) {
    ERROR_COND(base < 0 || base >= 256,
               "Non-standard bases were encountered in either the input reference or BAM file(s)");
    return kBaseIndexMap[base];
  }
 private:

 public:
const static std::vector<char> kBases;
const static int kBaseIndexMap[256];
};

} // easehts
} // ncic

#endif
