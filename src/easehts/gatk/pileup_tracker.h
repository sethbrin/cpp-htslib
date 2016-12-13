//
// Created by zp on 12/8/16.
//

#ifndef EASEHTSLIB_PILEUP_TRACKER_H_
#define EASEHTSLIB_PILEUP_TRACKER_H_

#include "easehts/gatk/alignment_state_machine.h"

#include <htslib/sam.h>

namespace ncic {
namespace easehts {
namespace gatk {

class PileupTracker {
 public:
  int start;
  int end;
  AlignmentStateMachine state_machine;
  SAMBAMRecord* read;

  PileupTracker(SAMBAMRecord* read)
    : state_machine(read) {
      this->read = read;
      start = read->GetAlignmentStart();
      end = read->GetAlignmentEnd();
  }

  ~PileupTracker() {
    delete read;
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

} // gatk
} // easehts
} // ncic

#endif
