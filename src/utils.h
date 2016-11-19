//
// Created by zp on 8/24/16.
//

#ifndef EASEHTSLIB_UTILS_H
#define EASEHTSLIB_UTILS_H

#include <htslib/sam.h>
#include <memory>
#include <string>
#include <stdio.h>
#include <vector>
#include <sys/stat.h>
// store some common values

#include <assert.h>

# define ASSERT(expr) assert(expr)
# define ASSERT_MSG(expr, msg) assert((expr)&&(msg))

#define CHECK_NOTNULL(val) ASSERT_MSG(val != NULL, "'" #val "' Must be non NULL")

#define ERROR_COND(cond, msg) do { if ((cond)) { fprintf(stderr, "[E::%s] %s \n", __func__, msg);  ::abort();} } while (0)
#define WARN_COND(cond, msg) if ((cond)) fprintf(stderr, "[W::%s] %s \n", __func__, msg)
#define ERROR(msg) do { fprintf(stderr, "[E::%s] %s \n", __func__, msg);  ::abort(); } while (0)
#define WARN(msg) fprintf(stderr, "[W::%s] %s \n", __func__, msg)

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

template<typename ... Args>
std::string StringFormat( const std::string& format, Args ... args )
{
  size_t size = snprintf(nullptr, 0, format.c_str(), args ...) + 1; // Extra space for '\0'
  std::unique_ptr<char[]> buf( new char[ size ]);
  snprintf(buf.get(), size, format.c_str(), args ...);
  return std::string(buf.get(), buf.get() + size - 1); // We don't want the '\0' inside
}

template<typename ... Args>
const char* StringFormatCStr( const std::string& format, Args ... args )
{
  size_t size = snprintf(nullptr, 0, format.c_str(), args ...) + 1; // Extra space for '\0'
  std::unique_ptr<char[]> buf( new char[ size ]);
  snprintf(buf.get(), size, format.c_str(), args ...);
  return std::string(buf.get(), buf.get() + size - 1).c_str(); // We don't want the '\0' inside
}

inline bool EndsWith(const std::string& value, const std::string& ending)
{
  if (ending.size() > value.size()) return false;
  return std::equal(ending.rbegin(), ending.rend(), value.rbegin());
}

inline bool FileExists(const std::string filename) {
  struct stat buffer;
  return (stat(filename.c_str(), &buffer) == 0);
}

} // utils
} // easehts
} // ncic
#endif //EASEHTSLIB_UTILS_H
