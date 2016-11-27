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

  double GetLikelihood(easehts::DiploidGenotype g) const {
    return GetLikelihoods()[g.GetIndex()];
  }


  easehts::DiploidSNPGenotypeLikelihoods CalculateGenotypeLikelihoods(
      char observed_base1, char quality_score1,
      char observed_base2, char quality_score2) {
    std::vector<double> log10_four_base_likelihoods =
      ComputeLog10Likelihoods(observed_base1, quality_score1,
                              observed_base2, quality_score2);

    VariableAllelicRatioGenotypeLikelihoods gl = *this;
    gl.SetToZero();

    for (const auto& g : easehts::DiploidGenotype::kGenotypes) {
      double f_base1 = 0;
      double f_base2 = 0;

      if (g.GetBase1() == ref_ || g.GetBase2() == ref_) {
        // if it is ref/ref use half
        if (g.GetBase1() == g.GetBase2()) {
          f_base1 = log_half_;
          f_base2 = log_half_;
        } else {
          // if one base is reference, use f
          f_base1 = (g.GetBase1() == ref_) ? log_one_minus_F_ : log_F_;
          f_base2 = (g.GetBase2() == ref_) ? log_one_minus_F_ : log_F_;
        }
      } else {
        f_base1 = log_half_;
        f_base2 = log_half_;
      }
      double p_base = 0.0;
      p_base += std::pow(10,
                         log10_four_base_likelihoods[easehts::BaseUtils::SimpleBaseToBaseIndex(g.GetBase1())]
                         - f_base1);
      p_base += std::pow(10,
                         log10_four_base_likelihoods[easehts::BaseUtils::SimpleBaseToBaseIndex(g.GetBase2())]
                         - f_base2);

      gl.log10_likelihoods_[g.GetIndex()] += std::log10(p_base);

    }
    return gl;
  }

  int Add(char observed_base1, char quality_score1) {
    char observed_base2 = 0;
    char quality_score2 = 0;

    // TODO current not just cache, so here not complete it
    // for more information, refer to the java code
    DiploidSNPGenotypeLikelihoods gl = CalculateGenotypeLikelihoods(
        observed_base1, quality_score1, observed_base2, quality_score2);

    const std::vector<double>& likelihoods = gl.GetLikelihoods();
    for (const auto& g : easehts::DiploidGenotype::kGenotypes) {
      log10_likelihoods_[g.GetIndex()] += likelihoods[g.GetIndex()];
    }
    return 1;
  }

  int Add(const easehts::PileupElement& element, bool ignore_bad_bases,
          bool cap_base_quals_at_mapping_qual, int min_base_qual) {
    char base = element.GetBase();
    char qual = QualToUse(element, ignore_bad_bases,
                          cap_base_quals_at_mapping_qual, min_base_qual);
    return qual == 0 ? 0 : Add(base, qual, 0, 0, 1);
  }


  /**
   *
   * @param observed_base1    first observed base
   * @param qual1       base qual of first observed base
   * @param observed_base2    second observed base
   * @param qual2       base qual of second observed base; can be 0, indicating no second base was observed for this fragment
   * @param n_obs        the number of times this quad of values was seen.  Generally 1, but reduced reads can have nObs > 1 for synthetic reads
   * @return 0 if the base is bad, 1 otherwise
   */

  int Add(
    char observed_base1, char qual1,
    char observed_base2, char qual2, int n_obs) {
  // TODO: current no cache
  DiploidSNPGenotypeLikelihoods gl =
    CalculateGenotypeLikelihoods(observed_base1, qual1,
                                 observed_base2, qual2);
  const std::vector<double>& likelihoods = gl.GetLikelihoods();
  for (const auto& g : easehts::DiploidGenotype::kGenotypes) {
    log10_likelihoods_[g.GetIndex()] += likelihoods[g.GetIndex()] * n_obs;
  }
  return 1;

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
