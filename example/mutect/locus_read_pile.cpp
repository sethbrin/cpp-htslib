#include "locus_read_pile.h"

#include <easehts/pileup.h>
#include <easehts/base_utils.h>
#include <easehts/diploid_genotype.h>

#include <vector>
#include <cmath>

namespace ncic {
namespace mutect {

const int LocusReadPile::kGapEventProximity = 5; // a 11bp window

LocusReadPile::LocusReadPile(SampleType sample_type, char ref_base, int min_quality_score,
                             int min_qsum_quality_score, bool allow_mapq0_for_qual_sum,
                             bool retain_overlap_mismatches, bool track_base_quality_socres)
  : sample_type_(sample_type),
    ref_base_(ref_base),
    min_quality_score_(min_quality_score),
    min_qsum_quality_score_(min_qsum_quality_score),
    allow_mapq0_for_qual_sum_(allow_mapq0_for_qual_sum),
    retain_overlap_mismatches_(retain_overlap_mismatches),
    track_base_quality_socres_(track_base_quality_socres),
    quality_sums_(track_base_quality_socres),
    deletions_count_(0),
    insertion_count_(0) {

}


void LocusReadPile::AddPileupElement(const easehts::ReadBackedRawPileup& read_backed_pileup) {
  for (size_t i=0; i < read_backed_pileup.Size(); i++) {
    pileup_.AddElement(read_backed_pileup[i]);
  }
}

void LocusReadPile::InitPileups() {
  easehts::ReadBackedPileup no_overlap_pileup;
  pileup_.GetOverlappingFragmentFilteredPileup(&no_overlap_pileup,
                                                       ref_base_, retain_overlap_mismatches_);

  no_overlap_pileup.GetPileupWithoutDeletions(&initial_pileup_);
  initial_pileup_.GetBaseFilteredPileup(min_quality_score_,
                                        &quality_score_filter_pileup_);
  quality_score_filter_pileup_.GetPileupWithoutMappingQualityZeroReads(&final_pileup_);

  easehts::ReadBackedPileup tmp_pileup;

  // GetPileupWithoutDeletions && GetBaseFilteredPileup
  auto pred = [this](easehts::PileupElement element)->bool {
    return !element.IsDeletion() &&
      (element.IsDeletion() || element.GetQual() >= this->min_quality_score_);
  };
  pileup_.GetPileupByFilter(&tmp_pileup, pred);

  tmp_pileup.GetPositiveStrandPileup(&final_pileup_positive_strand_);
  tmp_pileup.GetNegativeStrandPileup(&final_pileup_negative_strand_);

  for (size_t idx=0; idx<quality_score_filter_pileup_.Size(); idx++) {
    const easehts::PileupElement& p = quality_score_filter_pileup_[idx];
    if (p.GetMappingQuality() == 0 &&
        !allow_mapq0_for_qual_sum_) continue;
    if (p.GetQual() <= min_qsum_quality_score_) continue;
    if (p.GetQual() > min_qsum_quality_score_) {
      quality_sums_.IncrementSum(p.GetBase(), 1, p.GetQual());
    }
  }


  // Calculate how many are at this site and how many insertion
  // are within INSERTION_PROXIMITY bp
  for (size_t idx=0; idx<quality_score_filter_pileup_.Size(); idx++) {
    const easehts::PileupElement& p = quality_score_filter_pileup_[idx];
    if (p.GetBase() == easehts::PileupElement::kDeletionBase) {
      deletions_count_++;
    } else {
      // check for nearby events
      std::vector<easehts::CigarElement> cigars =
        easehts::SAMBAMRecord::ParseRawCigar(p.GetRead());
      int event_start = 0;
      for (const easehts::CigarElement& cigar : cigars) {
        if (cigar.GetOperator() == easehts::CigarElement::INSERTION &&
            std::abs(event_start - p.GetOffset()) < kGapEventProximity) {
          insertion_count_ ++;
          break;
        }

        if (cigar.GetOperator() == easehts::CigarElement::DELETION &&
            std::abs(event_start - p.GetOffset()) < kGapEventProximity) {
          deletions_count_ ++;
          break;
        }

        event_start += cigar.GetLength();
      }
    }
  }
}


double LocusReadPile::EstimateAlleleFraction(char ref, char alt) const {
  return LocusReadPile::EstimateAlleleFraction(final_pileup_, ref, alt);
}

double LocusReadPile::EstimateAlleleFraction(
    const easehts::ReadBackedPileup& read_backed_pileup,
    char ref, char alt) {
  std::vector<int> counts = read_backed_pileup.GetBaseCounts();
  int ref_count = counts[easehts::BaseUtils::SimpleBaseToBaseIndex(ref)];
  int alt_count = counts[easehts::BaseUtils::SimpleBaseToBaseIndex(alt)];

  int depth = ref_count + alt_count;

  return depth == 0 ? 0 : (static_cast<double>(alt_count)/ depth);
}

double LocusReadPile::CalculateLogLikelihood(
    const easehts::ReadBackedPileup& read_backed_pileup,
    char ref, char alt, double f) {
  double ll = 0;
  for (int i = 0; i < read_backed_pileup.Size(); i++) {
    const easehts::PileupElement& element = read_backed_pileup[i];

    char base = element.GetBase();
    char qual = element.GetQual();

    double e = std::pow(10, (qual / -10.0));

    if (base == ref) {
      ll += std::log10(f * e / 3 + (1 - f) * (1 - e));
    } else if (base == alt) {
      ll += std::log10(f * (1 - e) + (1 - f) * e / 3);
    } else {
      ll + std::log10(2 * e / 3);
    }
  }
  return ll;
}

std::array<double, 3> LocusReadPile::ExtractRefHemHom(
    const VariableAllelicRatioGenotypeLikelihoods& gl,
    char ref_allele, char alt_allele) {
  double ref = 0;
  double het = 0;
  double hom = 0;

  for (const auto& gt : easehts::DiploidGenotype::kGenotypes) {
    double likelihood = gl.GetLikelihood(gt);

    if (gt.GetBase1() == ref_allele &&
        gt.GetBase2() == ref_allele) {
      ref = likelihood;
    }

    if ((gt.GetBase1() == ref_allele &&
        gt.GetBase2() == alt_allele) ||
        (gt.GetBase1() == alt_allele &&
         gt.GetBase2() == ref_allele)) {
      het = likelihood;
    }

    if (gt.GetBase1() == alt_allele &&
        gt.GetBase2() == alt_allele) {
      hom = likelihood;
    }
  }
  return {ref, het, hom};
}


const easehts::DiploidGenotype& LocusReadPile::GetBestGenotype(
    const VariableAllelicRatioGenotypeLikelihoods& likelihoods) const {
  int idx = 0;
  int best_idx = 0;
  double best_likelihood = 0;
  for (const auto& gt : easehts::DiploidGenotype::kGenotypes) {
    double likelihood = likelihoods.GetLikelihood(gt);
    if (likelihood >= best_likelihood) {
      best_idx = idx;
      best_likelihood = likelihood;
    }
    idx ++;
  }
  return easehts::DiploidGenotype::kGenotypes[best_idx];
}

double LocusReadPile::GetRefVsAlt(
    const VariableAllelicRatioGenotypeLikelihoods& likelihoods,
    char ref, char alt_allele) {
  std::array<double, 3> ref_het_hom = ExtractRefHemHom(likelihoods, ref, alt_allele);

  return ref_het_hom[0] - LogAddSafe(ref_het_hom[1], ref_het_hom[2]);
}

double LocusReadPile::LogAddSafe(double a, double b) {
  double max_one = std::max(a, b);
  double min_one = std::min(a, b);

  return max_one + std::log(1 + std::pow(10, min_one - max_one));
}

double LocusReadPile::GetHetVsRef(
    const VariableAllelicRatioGenotypeLikelihoods& likelihoods,
    char ref, char alt_allele) {
  std::array<double, 3> ref_het_hom = ExtractRefHemHom(likelihoods, ref, alt_allele);

  return ref_het_hom[1] - ref_het_hom[0];
}

double LocusReadPile::GetAltVsRef(
    const VariableAllelicRatioGenotypeLikelihoods& likelihoods,
    char ref, char alt_allele) {
  std::array<double, 3> ref_het_hom = ExtractRefHemHom(likelihoods, ref, alt_allele);

  return LogAddSafe(ref_het_hom[1], ref_het_hom[2]) - ref_het_hom[0];
}

VariableAllelicRatioGenotypeLikelihoods LocusReadPile::CalculateLikelihoods(
    const easehts::ReadBackedPileup& pileup) const {
  return CalculateLikelihoods(0.5, pileup);
}

VariableAllelicRatioGenotypeLikelihoods LocusReadPile::CalculateLikelihoods(
    double alpha, const easehts::ReadBackedPileup& pileup) const {

  VariableAllelicRatioGenotypeLikelihoods likelihoods(ref_base_, alpha);
  // we have to do this rather than pass in the ReadBackedPileup because that call
  // attempts to make Fragments out of these, which doesn't work if you have
  // a single pileup with multiple samples (as we do in the simulation)
  int size = pileup.Size();
  for (int i=0; i<size; i++) {
    likelihoods.Add(pileup[i], false, false, min_quality_score_);
  }
  return likelihoods;
}




} // mutect
} // ncic
