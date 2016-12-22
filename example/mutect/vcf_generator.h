//
// Created by zp on 12/17/16.
//

#ifndef MUTECT_VCF_GENERATOR_H_
#define MUTECT_VCF_GENERATOR_H_

#include "candidate_mutation.h"

#include <easehts/noncopyable.h>
#include <easehts/vcf.h>
#include <easehts/sam_sequence_dictionary.h>

#include <string>

namespace ncic {
namespace mutect {

class VCFGenerator : public easehts::NonCopyable {
 public:
  VCFGenerator(const std::string& filename)
    : writer_(filename)
  {}

  void WriteHeader(const easehts::SAMSequenceDictionary& dict,
                   const std::string assemby);

  void Add(const CandidateMutation& candidate);

  static const std::string kTumorName;
  static const std::string kNormalName;

 private:
  easehts::VCFWriter writer_;
};

} // mutect
} // ncic

#endif
