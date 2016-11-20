//
// Created by zp on 9/3/16.
//

#include <string>
#include <vector>

#include "genome_loc.h"
#include "noncopyable.h"

namespace ncic {

namespace easehts {

class Interval {
 public:
  Interval(std::string sequence, int start, int end, bool negative, std::string name)
    : sequence_(sequence),
    start_(start),
    end_(end) {}

  Interval()
 private:
  std::string sequence_;
  int start_;
  int end_;
  bool negative_strand_;
  std::string name_;
};

class IntervalUtils : public NonCopyable {
  //std::vector<GenomeLoc> parseIntervalFile()
};

} // easehts
} // ncic
