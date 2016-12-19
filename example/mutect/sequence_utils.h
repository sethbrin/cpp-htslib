//
// Created by zp on 11/27/16.
//
#ifndef MUTECT_SEQUCENCE_UTILS_H_
#define MUTECT_SEQUCENCE_UTILS_H_

#include <easehts/base_utils.h>
#include <easehts/genome_loc.h>
#include <easehts/noncopyable.h>
#include <easehts/sam_bam_record.h>
#include <easehts/gatk/pileup.h>

#include <array>
#include <cmath>

namespace ncic {
namespace mutect {

class SequenceUtils : public easehts::NonCopyable {
 public:

  static bool IsReadHeavilySoftClipped(easehts::SAMBAMRecord* b,
                                       float threshold) {
    int total = 0;
    int clipped = 0;
    const std::vector<easehts::CigarElement>& cigars =
      b->GetCigar();

    for (const auto& ce : cigars) {
      total += ce.GetLength();
      if (ce.GetOperator() == easehts::CigarElement::SOFT_CLIP) {
        clipped += ce.GetLength();
      }
    }

    float rate = static_cast<float>(clipped) / static_cast<float>(total);

    return rate >= threshold;
  }

  static std::array<int, 4> GetStrandContingencyTable(
      const easehts::gatk::ReadBackedPileup& forward_pileup,
      const easehts::gatk::ReadBackedPileup& reverse_pileup,
      char ref, char alt) {
    // Construct a 2x2 contingency table of
    //            forward     reverse
    //      REF    a       b
    //      MUT    c       d
    //
    // and return an array of {a,b,c,d}
    int ref_idx = easehts::BaseUtils::SimpleBaseToBaseIndex(ref);
    int alt_idx = easehts::BaseUtils::SimpleBaseToBaseIndex(alt);

    std::array<int, 4> forward_base_counts =
      forward_pileup.GetBaseCounts();
    std::array<int, 4> reverse_base_counts =
      reverse_pileup.GetBaseCounts();

    return {forward_base_counts[ref_idx], reverse_base_counts[ref_idx],
    forward_base_counts[alt_idx], reverse_base_counts[alt_idx]};
  }

  static std::vector<int> GetOffsetsInRead(
      const easehts::gatk::ReadBackedPileup& pileup,
      const easehts::GenomeLoc& location,
      bool use_forward_offsets) {
    std::vector<int> positions;
    positions.reserve(pileup.Size());

    for (int idx=0; idx<pileup.Size(); idx++) {
      const easehts::gatk::PileupElement* p = pileup[idx];
      int position = location.GetStart();
      if (use_forward_offsets) {
        position -= p->GetRead()->GetAlignmentStart();
      } else {
        position -= p->GetRead()->GetAlignmentEnd();
      }
      positions.push_back(std::abs(position));
    }
    return positions;
  }

  static std::vector<int> GetForwardOffsetsInRead(
      const easehts::gatk::ReadBackedPileup& pileup,
      const easehts::GenomeLoc& location) {
    return GetOffsetsInRead(pileup, location, true);
  }

  static std::vector<int> GetReverseOffsetsInRead(
      const easehts::gatk::ReadBackedPileup& pileup,
      const easehts::GenomeLoc& location) {
    return GetOffsetsInRead(pileup, location, false);
  }

};

} // mutect
} // ncic

#endif
