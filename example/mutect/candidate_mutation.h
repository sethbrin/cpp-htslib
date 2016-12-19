//
// Created by zp on 11/20/16.
//

#ifndef MUTECT_CANDIDATE_MUTATION_H_
#define MUTECT_CANDIDATE_MUTATION_H_

#include <easehts/genome_loc.h>
#include <easehts/vcf.h>
#include <easehts/diploid_genotype.h>

#include <vector>
#include <array>
#include <cfloat>

namespace ncic {
namespace mutect {

class CandidateMutation {
 public:
  CandidateMutation(const easehts::GenomeLoc& loc, char ref){
    location = loc;
    ref_allele = ref;
  }

  // NOTE not implemented
  int GetCountOfNormalsObservedIn() const {
    return 0;
  }

  void AddRejectionReason(std::string reason) {
    rejected = true;
    rejection_reasons.push_back(reason);
  }

  bool IsSeenInPanelOfNormals() const {
    // TODO add VC
    return false;
  }

  bool IsGermlineAtRisk() const {
    // TODO
    return dbsnp_site && !cosmic_site;
  }

  bool IsRejected() const {
    return rejected;
  }

  // TODO
  int GetCountsOfNormalObservedIn() const {
    return 0;
  }

  easehts::GenomeLoc location;
  std::string sequence_context;
  char ref_allele;
  bool dbsnp_site = false;
  bool cosmic_site = false;

  // NOTE current not support
  easehts::VariantContext* panel_of_normals_VC;
  easehts::VariantContext* dbsnp_VC;

  bool covered = false;
  double power = 0;
  double tumor_power = 0;
  double normal_power = 0;
  double normal_power_with_snp_prior = 0;
  double normal_power_no_snp_prior = 0;

  char alt_allele = 'N';
  std::string tumor_sample_name = "TUMOR";
  std::string normal_sample_name = "NORMAL";

  double contamination_fraction = 0;

  double contaminant_lod = 0;
  int score = 0;

  int tumor_q20_count = 0;
  int normal_q20_count = 0;

  int total_reads = 0;
  int map_q0_reads = 0;
  int initial_tumor_ref_counts = 0;
  int initial_tumor_alt_counts = 0;
  int initial_tumor_ref_quality_sum = 0;
  int initial_tumor_alt_quality_sum = 0;
  int initial_tumor_non_ref_quality_sum = 0;
  int initial_tumor_read_depth = 0;
  int initial_normal_ref_counts = 0;
  int initial_normal_alt_counts = 0;
  int initial_normal_ref_quality_sum = 0;
  int initial_normal_alt_quality_sum = 0;
  int tumor_ref_max_map_q = 0;
  int tumor_alt_max_map_q = 0;
  int initial_normal_read_depth = 0;
  easehts::DiploidGenotype initial_normal_best_genotype;

  double initial_tumor_lod = 0;
  double initial_normal_lod = 0;

  double tumor_F = 0;
  double tumor_lodF_star = 0;
  double tumor_lodF_star_forward = 0;
  double tumor_lodF_star_reverse = 0;

  double normal_F = 0;

  double power_to_detect_positive_strand_artifact = 0;
  double power_to_detect_negative_strand_artifact = 0;

  std::array<int, 4> strand_contingency_table;

  std::vector<int> tumor_alt_forward_offsets_in_read;
  std::vector<int> tumor_alt_reverse_offsets_in_read;

  double tumor_forward_offsets_in_read_median = kDoubleUnintialized;
  double tumor_forward_offsets_in_read_mad = kDoubleUnintialized;
  double tumor_reverse_offsets_in_read_median = kDoubleUnintialized;
  double tumor_reverse_offsets_in_read_mad = kDoubleUnintialized;

  int tumor_insertion_count = 0;
  int tumor_deletion_count = 0;

  std::vector<int> tumor_ref_quality_scores;
  std::vector<int> tumor_alt_quality_scores;
  std::vector<int> normal_ref_quality_scores;
  std::vector<int> normal_alt_quality_scores;

  std::vector<std::string> rejection_reasons;
  bool rejected = false;

  // use the double max value
  const static double kDoubleUnintialized;

};


} // mutect
} // ncic

#endif
