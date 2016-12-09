//
// Created by zp on 12/8/16.
//

#ifndef EASEHTSLIB_PILEUP_TRACKER_H_
#define EASEHTSLIB_PILEUP_TRACKER_H_

#include "alignment_state_machine.h"

#include <htslib/sam.h>

namespace ncic {
namespace easehts {

class PileupTracker {
 public:
  int start;
  int end;
  AlignmentStateMachine state_machine;
  bam1_t* read;

  PileupTracker(bam1_t* read)
    : state_machine(read) {
      this->read = read;
      start = SAMBAMRecord::GetAlignmentStart(read);
      end = SAMBAMRecord::GetAlignmentEnd(read);
  }

  ~PileupTracker() {
    bam_destroy1(read);
  }

  bool StepForwardOnGenome() {
    if (state_machine.StepForwardOnGenome() == 0) {
      return false;
    } else {
      return true;
    }
  }

  bool IsBeforeEnd(int cur) {
    return cur <= end;
  }

  bool IsAfterStart(int cur) {
    return cur >= start;
  }

  int GetGenomeOffset() {
    return state_machine.GetGenomeOffset();
  }

};

} // easehts
} // ncic

#endif
