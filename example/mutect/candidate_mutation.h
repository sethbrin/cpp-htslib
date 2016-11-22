//
// Created by zp on 11/20/16.
//

#ifndef MUTECT_CANDIDATE_MUTATION_H_
#define MUTECT_CANDIDATE_MUTATION_H_

#include <easehts/genome_loc.h>

namespace ncic {
namespace mutect {

class CandidateMutation {
 public:
  CandidateMutation(const easehts::GenomeLoc& location, char ref_allele){
    location_ = location;
    ref_allele_ = ref_allele;
  }

  easehts::GenomeLoc location_;
  char ref_allele_;
};

} // mutect
} // ncic

#endif
