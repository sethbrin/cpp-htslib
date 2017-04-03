//
// Created by zp on 11/27/16.
//
#ifndef MUTECT_CGA_ALIGNMENT_UTILS_H_
#define MUTECT_CGA_ALIGNMENT_UTILS_H_

#include <easehts/noncopyable.h>
#include <easehts/gatk/pileup.h>
#include <easehts/sam_bam_record.h>
#include <easehts/reference_sequence.h>

#include <cmath>
#include <string>

namespace ncic {
namespace mutect {

class CGAAlignmentUtils : public easehts::NonCopyable {
 public:

  /** Returns the number of mismatches in the pileup element within the given reference context.
   *
   * @param reference the reference
   * @param p       the pileup element
   * @param location     current location
   * @param ignore_target_site     if true, ignore mismatches at the target locus (i.e. the center of the window)
   * @param quality_sum_instead_of_mismatch_count if true, return the quality score sum of the mismatches rather than the count
   * @return the number of mismatches
   */
  static int MismatchesInRefWindow(
      const easehts::ReferenceSequence& ref_bases,
      const easehts::gatk::PileupElement* p,
      const easehts::GenomeLoc& location,
      const easehts::GenomeLoc& window,
      bool ignore_target_site,
      bool quality_sum_instead_of_mismatch_count=false) {
    int sum = 0;

    int window_start = window.GetStart();
    int window_stop = window.GetStop();
    easehts::SAMBAMRecord* read = p->GetRead();
    std::string read_bases = read->GetSequence();
    uint8_t* read_qualities = read->GetRawQuality();
    const std::vector<easehts::CigarElement>& cigars =
      read->GetCigar();

    int read_index = 0;
    int current_pos = read->GetAlignmentStart();
    int ref_index = std::max(0, current_pos - window_start);

    for (const auto& ce : cigars) {
      int cigar_element_length = ce.GetLength();
      switch (ce.GetOperator()) {
        case easehts::CigarElement::MATCH:
          for (int j = 0; j < cigar_element_length; j++, read_index++, current_pos++) {
            // are we pass the ref window?
            if (current_pos > window_stop) break;

            // are we before the ref window?
            if (current_pos < window_start) continue;

            char ref_char = ref_bases[ref_index++];

            // do we need to skip the target site?
            if (ignore_target_site && location.GetStart() == current_pos) continue;

            char read_char = read_bases[read_index];
            if (ref_char != read_char) {
              sum += (quality_sum_instead_of_mismatch_count) ?
                read_qualities[read_index] : 1;
            }
          }
          break;
        case easehts::CigarElement::INSERTION:
        case easehts::CigarElement::SOFT_CLIP:
          read_index += cigar_element_length;
          break;
        case easehts::CigarElement::DELETION:
        case easehts::CigarElement::SKIPPED_REGION:
          current_pos += cigar_element_length;
          if (current_pos > window_start) {
            ref_index += std::min(cigar_element_length, current_pos - window_start);
          }
          break;
        default:break;
      }
    }
    return sum;
  }
};

} // mutect
} // ncic

#endif
