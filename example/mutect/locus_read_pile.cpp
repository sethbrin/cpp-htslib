#include "locus_read_pile.h"

#include <easehts/pileup.h>
#include <vector>
#include <cmath>

namespace ncic {
namespace mutect {

const int LocusReadPile::kGapEventProximity = 5; // a 11bp window

LocusReadPile::LocusReadPile(SampleType sample_type, char ref_base, int min_quality_score,
                             int min_qsum_quality_score, bool allow_mapq0_for_qual_sum,
                             bool retain_overlap_mismatches, bool track_base_quality_socres)
  : sample_type_(sample_type),
    ref_base_(ref_base),
    min_quality_score_(min_quality_score),
    min_qsum_quality_score_(min_qsum_quality_score),
    allow_mapq0_for_qual_sum_(allow_mapq0_for_qual_sum),
    retain_overlap_mismatches_(retain_overlap_mismatches),
    track_base_quality_socres_(track_base_quality_socres),
    quality_sums_(track_base_quality_socres),
    deletions_count_(0),
    insertion_count_(0) {

}


void LocusReadPile::AddPileupElement(const easehts::ReadBackedRawPileup& read_backed_pileup) {
  for (size_t i=0; i < read_backed_pileup.Size(); i++) {
    pileup_.AddElement(read_backed_pileup[i]);
  }
}

void LocusReadPile::InitPileups() {
  easehts::ReadBackedPileup no_overlap_pileup;
  pileup_.GetOverlappingFragmentFilteredPileup(&no_overlap_pileup,
                                                       ref_base_, retain_overlap_mismatches_);

  no_overlap_pileup.GetPileupWithoutDeletions(&initial_pileup_);
  initial_pileup_.GetBaseFilteredPileup(min_quality_score_,
                                        &quality_score_filter_pileup_);
  quality_score_filter_pileup_.GetPileupWithoutMappingQualityZeroReads(&final_pileup_);

  easehts::ReadBackedPileup tmp_pileup;

  // GetPileupWithoutDeletions && GetBaseFilteredPileup
  auto pred = [this](easehts::PileupElement element)->bool {
    return !element.IsDeletion() &&
      (element.IsDeletion() || element.GetQual() >= this->min_quality_score_);
  };
  pileup_.GetPileupByFilter(&tmp_pileup, pred);

  tmp_pileup.GetPositiveStrandPileup(&final_pileup_positive_strand_);
  tmp_pileup.GetNegativeStrandPileup(&final_pileup_negative_strand_);

  for (size_t idx=0; idx<quality_score_filter_pileup_.Size(); idx++) {
    const easehts::PileupElement& p = quality_score_filter_pileup_[idx];
    if (p.GetMappingQuality() == 0 &&
        !allow_mapq0_for_qual_sum_) continue;
    if (p.GetQual() <= min_qsum_quality_score_) continue;
    if (p.GetQual() > min_qsum_quality_score_) {
      quality_sums_.IncrementSum(p.GetBase(), 1, p.GetQual());
    }
  }


  // Calculate how many are at this site and how many insertion
  // are within INSERTION_PROXIMITY bp
  for (size_t idx=0; idx<quality_score_filter_pileup_.Size(); idx++) {
    const easehts::PileupElement& p = quality_score_filter_pileup_[idx];
    if (p.GetBase() == easehts::PileupElement::kDeletionBase) {
      deletions_count_++;
    } else {
      // check for nearby events
      std::vector<easehts::CigarElement> cigars =
        easehts::SAMBAMRecord::ParseRawCigar(p.GetRead());
      int event_start = 0;
      for (const easehts::CigarElement& cigar : cigars) {
        if (cigar.GetOperator() == easehts::CigarElement::INSERTION &&
            std::abs(event_start - p.GetOffset()) < kGapEventProximity) {
          insertion_count_ ++;
          break;
        }

        if (cigar.GetOperator() == easehts::CigarElement::DELETION &&
            std::abs(event_start - p.GetOffset()) < kGapEventProximity) {
          deletions_count_ ++;
          break;
        }

        event_start += cigar.GetLength();
      }
    }
  }
}

} // mutect
} // ncic
