//
// Created by zp on 11/21/16.
//

#ifndef MUTECT_LOCUS_READ_PILE_H_
#define MUTECT_LOCUS_READ_PILE_H_

#include "quality_sums.h"

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
  double EstimateAlleleFraction(char ref, char alt) const;

  double CalculateAltVsRefLOD(const easehts::ReadBackedPileup& pileup,
                              char alt, double alternate, double reference) const {
    return CalculateAltVsRefLOD(pileup, ref_base_, alt, alternate, reference);
  }

  double CalculateAltVsRefLOD(char alt, double alternate, double reference) const {
    return CalculateAltVsRefLOD(final_pileup_, alt, alternate, reference);
  }

  double CalculateAltVsRefLOD(const easehts::ReadBackedPileup& pileup, char ref,
                              char alt, double alternate, double reference) const {
    double lod_alt = LocusReadPile::CalculateLogLikelihood(pileup, ref, alt, alternate);
    double lod_ref = LocusReadPile::CalculateLogLikelihood(pileup, ref, alt, reference);

    return lod_alt - lod_ref;
  }

 public:
  const static int kGapEventProximity;

  static double EstimateAlleleFraction(const easehts::ReadBackedPileup& read_backed_pileup,
                                char ref, char alt);
  static double CalculateLogLikelihood(const easehts::ReadBackedPileup& read_backed_pileup,
                                       char ref, char alt, double f);

 public:
  easehts::ReadBackedPileup pileup_;
  easehts::ReadBackedPileup initial_pileup_;
  easehts::ReadBackedPileup quality_score_filter_pileup_;
  easehts::ReadBackedPileup final_pileup_;
  easehts::ReadBackedPileup final_pileup_positive_strand_;
  easehts::ReadBackedPileup final_pileup_negative_strand_;

  QualitySums quality_sums_;

  SampleType sample_type_;
  char ref_base_;
  int min_quality_score_;
  int min_qsum_quality_score_;
  bool allow_mapq0_for_qual_sum_;
  bool retain_overlap_mismatches_;
  bool track_base_quality_socres_;
  int deletions_count_;
  int insertion_count_;
};

} // mutect
} // ncic

#endif

