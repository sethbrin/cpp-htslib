//
// Created by zp on 11/25/16.
//

#ifndef EASEHTSLIB_DIPLOID_SNP_GENOTYPE_LIKELIDOODS_H_
#define EASEHTSLIB_DIPLOID_SNP_GENOTYPE_LIKELIDOODS_H_


#include <cmath>
#include <vector>

namespace ncic {
namespace easehts {


/**
 * Stable, error checking version of the Bayesian genotyper.  Useful for calculating the likelihoods, priors,
 * and posteriors given a pile of bases and quality scores
 *
 * Suppose we have bases b1, b2, ..., bN with qualities scores q1, q2, ..., qN.  This object
 * calculates:
 *
 * P(G | D) = P(G) * P(D | G)
 *
 * where
 *
 * P(D | G) = sum_i log10 P(bi | G)
 *
 * and
 *
 * P(bi | G) = 1 - P(error | q1) if bi is in G
 *           = P(error | q1) / 3 if bi is not in G
 *
 * for homozygous genotypes and for heterozygous genotypes:
 *
 * P(bi | G) = 1 - P(error | q1) / 2 + P(error | q1) / 6 if bi is in G
 *           = P(error | q1) / 3 if bi is not in G
 *
 * for each of the 10 unique diploid genotypes AA, AC, AG, .., TT
 *
 * Everything is stored as arrays indexed by DiploidGenotype.ordinal() values in log10 space.
 *
 * The priors contain the relative probabilities of each genotype, and must be provided at object creation.
 * From then on, you can call any of the add() routines to update the likelihoods and posteriors in the above
 * model.
 */
class DiploidSNPGenotypeLikelihoods {
 public:

  /**
   * Create a new GenotypeLikelhoods object with given PCR error rate for each diploid genotype
   *
   * @param PCR_error_rate  the PCR error rate
   */
  DiploidSNPGenotypeLikelihoods(double PCR_error_rate) {
    log10_PCR_error_3_ = std::log10(PCR_error_rate) -
      DiploidSNPGenotypeLikelihoods::kLog10_3;
    log10_1_minus_PCR_error_ = std::log10(1.0 - PCR_error_rate);
  }

  const std::vector<double>& GetLikelihoods() const {
    return log10_likelihoods_;
  }

  void SetToZero() {
    log10_likelihoods_.clear();
  }

  static const double kDefaultPcrErrorRate;

 protected:
  DiploidSNPGenotypeLikelihoods CalculateGenotypeLikelihoods(
      char observed_base1, char quality_score1,
      char observed_base2, char quality_score2);

  /**
   * Updates likelihoods and posteriors to reflect an additional observation of observedBase with
   * qualityScore.
   *
   * @param observedBase1  the base observed on the 1st read of the fragment
   * @param qualityScore1  the qual of the base on the 1st read of the fragment, or zero if NA
   * @param observedBase2  the base observed on the 2nd read of the fragment
   * @param qualityScore2  the qual of the base on the 2nd read of the fragment, or zero if NA
   * @return likelihoods for this observation or null if the base was not considered good enough
   * to add to the likelihoods (Q0 or 'N', for example)
   */
  std::vector<double> ComputeLog10Likelihoods(
      char observed_base1, char quality_score1,
      char observed_base2, char quality_score2);

  static double Log10PofObservingBaseGivenChromosome(
      char observed_base,char chrom_base, char qual);

  //
  // The fundamental data arrays associated with a Genotype Likelihoods object
  //
  std::vector<double> log10_likelihoods_;

  // TODO: don't calulate this each time through
  double log10_PCR_error_3_;
  double log10_1_minus_PCR_error_;

 protected:
  static const int kFixedPloidy;
  static const int kMaxPloidy;
  static const double kPloidyAdjustment;
  static const double kLog10_3;

 private:
};

} // easehts
} // ncic

#endif
