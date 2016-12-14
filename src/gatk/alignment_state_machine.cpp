#include "easehts/gatk/alignment_state_machine.h"
#include "easehts/gatk/pileup.h"

namespace ncic {
namespace easehts {
namespace gatk {

char AlignmentStateMachine::StepForwardOnGenome() {
  // loop util we either find a cigar element step that moves us one base on
  // the genome, or we run out of cigar elements
  const std::vector<CigarElement>& cigars = read_->GetCigar();
  while (true) {
    // we enter this method with read_offset_ == index of the last
    // processed base ob the read(-1 if we did not process a single base yet);
    // this can be last maching base, or last base of an insertion
    if (current_cigar_element_offset_ == -1 ||
        (offset_into_current_cigar_element_ + 1) >=
        cigars[current_cigar_element_offset_].GetLength()) {
      current_cigar_element_offset_++;

      if (current_cigar_element_offset_ < n_cigars_) {
        offset_into_current_cigar_element_ = -1;
        current_element_ = cigars[current_cigar_element_offset_];
        // next line: guards against cigar element of length 0; when new cigar
        // element is retrieved, we reenter in order to re-check
        // offset_into_current_cigar_element_ against current element's length
        continue;
      } else {
        if (current_cigar_element_offset_ != -1 &&
            current_cigar_element_offset_ < n_cigars_ &&
            cigars[current_cigar_element_offset_].GetOperator() ==
            CigarElement::DELETION) {
          ERROR(utils::StringFormatCStr(
                  "read ends with deletion. Cigar:%s. "
                  "Although the SAM spec technically permits such reads,"
                  "this is often indicative of malformed files. "
                  "If you are sure you want to use this file, "
                  "re-run your analysis with the extra option:"
                  "-rf BadCigar", read_->GetCigarString().c_str()));
        }

        // we're done, so set the offset of the cigar to 0 for cleanliness,
        // as well as the current element
        offset_into_current_cigar_element_ = 0;
        current_cigar_element_offset_ = -1;
        read_offset_ = read_->GetSequenceLength();

        // Reads the contain indels model the genome_offset_ as the following
        // base in the reference. Because we fall into this elese block only
        // when indels end the read, increment genome_offset_ such
        // that the current offset of the read is the next ref base after the
        // end of the indel. This position will model a point on the
        // referece somewhere after the end of the read.

        // extend events need that. Logically, it's legal to advance the
        // genomic offset here:
        genome_offset_ ++;

        // we do step forward on the ref, and by returning 0 we also
        // indicate that we are past the read end.
        return 0;

      }
    }

    offset_into_current_cigar_element_++;
    bool done = false;
    switch (current_element_.GetOperator()) {
      case CigarElement::HARD_CLIP:
      case CigarElement::PADDING:
        offset_into_current_cigar_element_ = current_element_.GetLength();
        break;
      case CigarElement::INSERTION:
      case CigarElement::SOFT_CLIP:
        offset_into_current_cigar_element_ = current_element_.GetLength();
        read_offset_ += current_element_.GetLength();
        break;
      case CigarElement::DELETION:
        ERROR_COND(read_offset_ < 0,
                   utils::StringFormatCStr(
                       "read starts with deletion. Cigar: %s. "
                       "Although the SAM spec technically permits such reads,"
                       " this is often indicative of malformed files. "
                       "If you are sure you want to use this file,"
                       "re-run your analysis with the extra option: "
                       "-rf BadCigar", read_->GetCigarString().c_str()));
        // should be the same as N case
        genome_offset_++;
        done = true;
        break;
      case CigarElement::SKIPPED_REGION:
        genome_offset_ ++;
        done = true;
        break;
      case CigarElement::MATCH:
      case CigarElement::EQUAL:
      case CigarElement::MISMATCH:
        read_offset_++;
        genome_offset_++;
        done = true;
        break;
      default:
        ERROR(utils::StringFormatCStr(
                "Case statement didn't deal with cigar op:%c",
                current_element_.GetOperator()));
    }

    if (done) return current_element_.GetOperator();
  }
}

bool AlignmentStateMachine::MakePileupElement(PileupElement** element) {
  if (IsLeftEdge() || IsRightEdge()) {
    WARN("Cannot make a pileup element from an edge alignment state");
    return false;
  }
  if (current_cigar_element_offset_ == -1 ||
      current_cigar_element_offset_ >= n_cigars_) {
    return false;
  }
  const std::vector<CigarElement>& cigars = read_->GetCigar();
  // N's are never added to any pileup
  // LocusIteratorByState.java
  if (cigars[current_cigar_element_offset_].GetOperator()
      == CigarElement::SKIPPED_REGION) {
    return false;
  }
  uint8_t qual = read_->GetRawQuality()[read_offset_];
  //uint8_t base = utils::SAMUtils::GetSequenceBaseChar(
  //      read_->GetRawSequence(), read_offset_);
  uint8_t base = read_->GetSequence()[read_offset_];

  bool is_del = cigars[current_cigar_element_offset_].GetOperator()
    == CigarElement::DELETION;
  *element = new PileupElement(read_, qual, base, is_del, read_offset_);

  return true;
}

} // gatk
} // easehts
} // ncic
