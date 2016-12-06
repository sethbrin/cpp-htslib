//
// Created by zp on 12/3/16.
//

#ifndef EASEHTSLIB_ALIGNMENT_STATE_MACHINE_H_
#define EASEHTSLIB_ALIGNMENT_STATE_MACHINE_H_

#include "noncopyable.h"
#include "sam_bam_record.h"
#include "sam_bam_reader.h"
#include "utils.h"

#include <vector>

namespace ncic {
namespace easehts {

class PileupElement;
class AlignmentStateMachine : public NonCopyable {
 public:
  AlignmentStateMachine(bam1_t* read) {
    read_ = read;
    cigars_ = SAMBAMRecord::ParseRawCigar(read_);
    n_cigars_ = cigars_.size();
    InitializeAsLeftEdge();
  }

  /**
   * Get the read we are aligning to the genome
   * @return a NONULL record
   */
  bam1_t* GetRawRead() const {
    return read_;
  }

  /**
   * Is this the left edge state? I.e., one that is before of after the current
   * read?
   *
   * @return true if this state is an edge state, false otherwise
   */
  bool IsLeftEdge() const {
    return read_offset_ == -1;
  }

  /**
   * Are we on the right edge? I.e., is the current state off the right of the
   * alignment?
   * @return true if off the right edge, false if otherwise
   */
  bool IsRightEdge() const {
    return read_offset_ == SAMBAMRecord::GetSequenceLength(read_);
  }

  /**
   * What is our current offset in the read's bases that aligns us with the
   * reference genome?
   *
   * @return the current read offset position? if an edge will be == -1
   */
  int GetReadOffset() const {
    return read_offset_;
  }

  /**
   * Get the offset of the current cigar element among all cigar elements in the read
   *
   * Suppose our read's cigar is 1M2D3M, and we're at the first 1M.  This would
   * return 0.  Stepping forward puts us in the 2D, so our offset is 1.  Another
   * step forward would result in a 1 again (we're in the second position of the 2D).
   * Finally, one more step forward brings us to 2 (for the 3M element)
   *
   * @return the offset of the current cigar element in the reads's cigar.  Will return -1 for
   * when the state is on the left edge, and be == the number of cigar elements in the
   * read when we're past the last position on the genome
   */
  int GetCurrentCigarElementOffset() const {
    return current_cigar_element_offset_;
  }

  /**
   * Get the offset of the current state into the current cigar element
   *
   * That is, suppose we have a read with cigar 2M3D4M, and we're right at
   * the second M position.  offsetIntoCurrentCigarElement would be 1, as
   * it's two elements into the 2M cigar.  Now stepping forward we'd be
   * in cigar element 3D, and our offsetIntoCurrentCigarElement would be 0.
   *
   * @return the offset (from 0) of the current state in the current cigar element.
   *  Will be 0 on the right edge, and -1 on the left.
   */
  int GetOffsetIntoCurrentCigarElement() {
    return offset_into_current_cigar_element_;
  }

  /**
   * Step the state machine forward one unit
   *
   * Takes the current state of this machine, and advances the state until the next on-genome
   * cigar element (M, X, =, D) is encountered, at which point this function returns with the
   * cigar operator of the current element.
   *
   * Assumes that the AlignmentStateMachine is in the left edge state at the start, so that
   * stepForwardOnGenome() can be called to move the machine to the first alignment position.  That
   * is, the normal use of this code is:
   *
   * AlignmentStateMachine machine = new AlignmentStateMachine(read)
   * machine.stepForwardOnGenome()
   * // now the machine is at the first position on the genome
   *
   * When stepForwardOnGenome() advances off the right edge of the read, the state machine is
   * left in a state such that isRightEdge() returns true and returns null, indicating the
   * the machine cannot advance further.  The machine may explode, though this is not contracted,
   * if stepForwardOnGenome() is called after a previous call returned null.
   *
   * @return the operator of the cigar element that machine stopped at, 0 if we advanced off the end of the read
   */
  char StepForwardOnGenome();

  /**
   * Create a new PileupElement based on the current state of this element
   *
   * Must not be a left or right edge
   *
   * @return a pileup element
   */
  bool MakePileupElement(PileupElement* element);

 private:
  /**
   * Initialize the state variables to put this machine one bp before the
   *start of the alignment, so that a call to stepForwardOnGenome()
   * will advance us to the first proper location
   */
  void InitializeAsLeftEdge() {
    read_offset_ = -1;
    offset_into_current_cigar_element_ = -1;
    genome_offset_ = -1;
  }

  bam1_t* read_;
  std::vector<CigarElement> cigars_;
  int n_cigars_;
  int current_cigar_element_offset_ = -1;
  CigarElement current_element_;

  /**
   * How far are we offset from start of the read bases?
   */
  int read_offset_;

  /**
   * How far are we offset from the alignment start on the genome?
   */
  int genome_offset_;

  /**
   * How far are we into our CigarElement
   */
  int offset_into_current_cigar_element_;

};

} // easehts
} // ncic

#endif
