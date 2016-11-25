#include "easehts/diploid_SNP_genotype_likelihoods.h"
#include "easehts/base_utils.h"
#include "easehts/diploid_genotype.h"

#include <cmath>

namespace ncic {
namespace easehts {

const double DiploidSNPGenotypeLikelihoods::kDefaultPcrErrorRate = 1e-4;
const int DiploidSNPGenotypeLikelihoods::kFixedPloidy = 2;
const int DiploidSNPGenotypeLikelihoods::kMaxPloidy = 3;
const double DiploidSNPGenotypeLikelihoods::kPloidyAdjustment = std::log10(kFixedPloidy);
const double DiploidSNPGenotypeLikelihoods::kLog10_3 = std::log10(3.0);

DiploidSNPGenotypeLikelihoods DiploidSNPGenotypeLikelihoods::CalculateGenotypeLikelihoods(
    char observed_base1, char quality_score1,
    char observed_base2, char quality_score2) {
  DiploidSNPGenotypeLikelihoods gl = *this;
  gl.SetToZero();

  // we need to adjust for ploidy.  We take the raw p(obs | chrom) / ploidy,
  // which is -log10(ploidy) in log space
  for (const auto& g : DiploidGenotype::kGenotypes) {
    double p_base = 0.0;
    //p_base += std::pow(10, )
  }
}

std::vector<double> DiploidSNPGenotypeLikelihoods::ComputeLog10Likelihoods(
    char observed_base1, char quality_score1,
    char observed_base2, char quality_score2) {
  int size = BaseUtils::kBases.size();
  std::vector<double> log10_four_base_likelihoods(size, 0);

  for (int i = 0; i < size; i++) {
    char true_base = BaseUtils::kBases[i];
    double likelihood = 0.0;

    for (int j = 0; j < size; j++) {
      char fragment_base = BaseUtils::kBases[j];
      double log10_fragment_likelihood = (true_base == fragment_base)
        ? log10_1_minus_PCR_error_ : log10_PCR_error_3_;

      if (quality_score1 != 0) {
        log10_fragment_likelihood += Log10PofObservingBaseGivenChromosome(
            observed_base1, fragment_base, quality_score1);
      } else {
        log10_fragment_likelihood += Log10PofObservingBaseGivenChromosome(
            observed_base2, fragment_base, quality_score2);
      }

      likelihood += std::pow(10, log10_fragment_likelihood);
    }
    log10_four_base_likelihoods[BaseUtils::kBaseIndexMap[true_base]] = std::log10(likelihood);
  }
}

double DiploidSNPGenotypeLikelihoods::Log10PofObservingBaseGivenChromosome(
    char observed_base,char chrom_base, char qual)
{
  double logP = 0;

  if (observed_base == chrom_base) {
    double e = std::pow(10, (qual / -10.0));
    logP = std::log10(1.0 - e);
  } else {
    logP = qual / -10.0 + (-kLog10_3);
  }

  return logP;
}


} // easehts
} // ncic
