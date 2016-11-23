//
// Created by zp on 11/23/16.
//

#ifndef MUTECT_QUALITY_SUM_H_
#define MUTECT_QUALITY_SUM_H_

#include <easehts/noncopyable.h>
#include <easehts/utils.h>

#include <vector>

namespace ncic {
namespace mutect {

class QualitySums : public easehts::NonCopyable {
 public:
  explicit QualitySums(bool enable_quality_score_tracking)
    : enable_quality_score_tracking_(enable_quality_score_tracking) {
      Reset();
  }

  int GetQualitySum(char base) {
    if (base == 'a' || base == 'A') return a_score_;
    if (base == 'c' || base == 'C') return c_score_;
    if (base == 'g' || base == 'G') return g_score_;
    if (base == 't' || base == 'T') return t_score_;
    ERROR(easehts::utils::StringFormatCStr("Unknown base:%c", base));
  }

  int GetCounts(char base) {
    if (base == 'a' || base == 'A') return a_count_;
    if (base == 'c' || base == 'C') return c_count_;
    if (base == 'g' || base == 'G') return g_count_;
    if (base == 't' || base == 'T') return t_count_;
    ERROR(easehts::utils::StringFormatCStr("Unknown base:%c", base));
  }

  void IncrementSum(char base, int count, int qual_sum) {
    if (base == 'a' || base == 'A') {
      a_score_ += qual_sum;
      a_count_ += count;
      if (enable_quality_score_tracking_) {
        a_quality_sum_.push_back(qual_sum);
      }
      return;
    }
    if (base == 'c' || base == 'C') {
      c_score_ += qual_sum;
      c_count_ += count;
      if (enable_quality_score_tracking_) {
        c_quality_sum_.push_back(qual_sum);
      }
      return;
    }
    if (base == 'g' || base == 'G') {
      g_score_ += qual_sum;
      g_count_ += count;
      if (enable_quality_score_tracking_) {
        g_quality_sum_.push_back(qual_sum);
      }
      return;
    }
    if (base == 't' || base == 'T') {
      t_score_ += qual_sum;
      t_count_ += count;
      if (enable_quality_score_tracking_) {
        t_quality_sum_.push_back(qual_sum);
      }
      return;
    }
  }

  int GetOtherQualities(char base) {
    int total = a_score_ + c_score_ + g_score_ + t_score_;
    if (base == 'a' || base == 'A') return total - a_score_;
    if (base == 'c' || base == 'C') return total - c_score_;
    if (base == 'g' || base == 'G') return total - g_score_;
    if (base == 't' || base == 'T') return total - t_score_;
    ERROR(easehts::utils::StringFormatCStr("Unknown base:%c", base));
  }

  const std::vector<int>& GetBaseQualityScors(char base) {
    if (base == 'a' || base == 'A') return a_quality_sum_;
    if (base == 'c' || base == 'C') return c_quality_sum_;
    if (base == 'g' || base == 'G') return g_quality_sum_;
    if (base == 't' || base == 'T') return t_quality_sum_;
    ERROR(easehts::utils::StringFormatCStr("Unknown base:%c", base));
  }

  void Reset() {
    a_score_ = c_score_ = g_score_ = t_score_ = 0;
    a_count_ = c_count_ = g_count_ = t_count_ = 0;
    a_quality_sum_.clear();
    c_quality_sum_.clear();
    g_quality_sum_.clear();
    t_quality_sum_.clear();
  }

 private:
  int a_score_;
  int c_score_;
  int g_score_;
  int t_score_;
  int a_count_;
  int c_count_;
  int g_count_;
  int t_count_;

  // used for tracking individual base quality score
  bool enable_quality_score_tracking_;
  std::vector<int> a_quality_sum_;
  std::vector<int> c_quality_sum_;
  std::vector<int> g_quality_sum_;
  std::vector<int> t_quality_sum_;

};

} // mutect
} // ncic

#endif
