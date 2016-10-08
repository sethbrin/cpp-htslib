//
// Created by zp on 8/24/16.
//

#ifndef EASEHTSLIB_UTILS_H
#define EASEHTSLIB_UTILS_H

#include <vector>
#include <string>
#include <htslib/sam.h>
// store some common values

#include <assert.h>
# define ASSERT(expr) assert(expr)
# define ASSERT_MSG(expr, msg) assert((expr)&&(msg))

#define CHECK_NOTNULL(val) ASSERT_MSG(val != NULL, "'" #val "' Must be non NULL")

namespace ncic {

namespace easehts {

namespace utils {

void tokenize(const std::string &s, char c, std::vector<std::string> *res);


class SAMUtils {
 public:
  /*! @function
   @abstract  Get a base on read
   @param  s  Query sequence returned by bam_get_seq()
   @param  i  The i-th position, 0-based
   @return    4-bit integer representing the base.
   */
  static uint8_t GetSequenceBase(const uint8_t* s, int idx) {
    return bam_seqi(s, idx);
  }

  /**
   * Get the char of the read
   */
  static char GetSequenceBaseChar(const uint8_t* s, int idx) {
    return "=ACMGRSVTWYHKDBN"[GetSequenceBase(s, idx)];
  }

};

} // util
} // easehts
} // ncic
#endif //EASEHTSLIB_UTILS_H
