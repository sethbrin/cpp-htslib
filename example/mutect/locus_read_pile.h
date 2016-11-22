//
// Created by zp on 11/21/16.
//

#ifndef MUTECT_LOCUS_READ_PILE_H_
#define MUTECT_LOCUS_READ_PILE_H_

#include <easehts/noncopyable.h>
#include <easehts/pileup.h>

#include <vector>

namespace ncic {
namespace mutect {

enum SampleType {
  TUMOR, NORMAL
};

class LocusReadPile : public easehts::NonCopyable {
 public:
  LocusReadPile(SampleType sample_type, char ref_base, int min_quality_score,
                int min_qsum_quality_score, bool allow_mapq0_for_qual_sum,
                bool retain_overlap_mismatches, bool track_base_quality_socres);

  LocusReadPile(SampleType sample_type, char ref_base, int min_quality_score,
                int min_qsum_quality_score, bool track_base_quality_socres)
    : LocusReadPile(sample_type, ref_base, min_quality_score,
                    min_qsum_quality_score, false, false, track_base_quality_socres) {}

  void AddPileupElement(const easehts::ReadBackedRawPileup& read_backed_pileup);

  // init the pileups
  void InitPileups();

  int Size() const {
    return pileup_.Size();
  }

 private:
  easehts::ReadBackedPileup pileup_;
  easehts::ReadBackedPileup initial_pileup_;
  easehts::ReadBackedPileup quality_score_filter_pileup_;
  easehts::ReadBackedPileup final_pileup_;
  easehts::ReadBackedPileup final_pileup_positive_strand_;
  easehts::ReadBackedPileup final_pileup_negative_strand_;

  SampleType sample_type_;
  char ref_base_;
  int min_quality_score_;
  int min_qsum_quality_score_;
  bool allow_mapq0_for_qual_sum_;
  bool retain_overlap_mismatches_;
  bool track_base_quality_socres_;
};

} // mutect
} // ncic

#endif

