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
  LocusReadPile(SampleType sample_type)
    : sample_type_(sample_type)
  {
  }

  void AddPileupElement(const easehts::ReadBackedPileup& read_backed_pileup);

  int Size() const {
    return elements_.size();
  }

  void Reset() {
    elements_.clear();
  }
 private:
  std::vector<easehts::PileupElement> elements_;
  SampleType sample_type_;
};

} // mutect
} // ncic

#endif

