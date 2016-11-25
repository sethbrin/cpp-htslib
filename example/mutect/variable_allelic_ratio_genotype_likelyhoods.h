//
// Created by zp on 11/25/16.
//

#ifndef MUTECT_VARIABLE_ALLELIC_RATIO_GENOTYPE_LKEELIHOODS_H_
#define MUTECT_VARIABLE_ALLELIC_RATIO_GENOTYPE_LKEELIHOODS_H_

#include "easehts/diploid_SNP_genotype_likelihoods.h"
#include "easehts/diploid_genotype.h"

#include <cmath>

namespace ncic {
namespace mutect {

class VariableAllelicRatioGenotypeLikelihoods : easehts::DiploidSNPGenotypeLikelihoods {
 public:
  VariableAllelicRatioGenotypeLikelihoods(char ref, double f)
    : easehts::DiploidSNPGenotypeLikelihoods(0) {
    ref_ = ref;

    log_F_ = std::log10(1 / f);
    log_one_minus_F_ = std::log10(1/(1-f));
    log_half_ = std::log10(1/0.5);
  }

  double GetLikelihood(DiploidGenotype g) {
    return GetLikelihoods()[g.GetIndex()];
  }

  double log_F_;
  double log_one_minus_F_;
  double log_half_;
  char ref_;
 private:
};

} // mutect
} // ncic


#endif
