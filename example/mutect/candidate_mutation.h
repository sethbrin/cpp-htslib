//
// Created by zp on 11/20/16.
//

#ifndef MUTECT_CANDIDATE_MUTATION_H_
#define MUTECT_CANDIDATE_MUTATION_H_

#include <easehts/genome_loc.h>
#include <easehts/variant_context.h>
#include <easehts/diploid_genotype.h>

#include <vector>

namespace ncic {
namespace mutect {

class CandidateMutation {
 public:
  CandidateMutation(const easehts::GenomeLoc& loc, char ref){
    location = loc;
    ref_allele = ref;
  }

  // NOTE not implemented
  int GetCountOfNormalsObservedIn() {
    return 0;
  }

  void AddRejectionReason(std::string reason) {
    rejected = true;
    rejection_reasons.push_back(reason);
  }

  easehts::GenomeLoc location;
  std::string sequence_context;
  char ref_allele;
  bool dbsnp_site = false;
  bool cosmic_site = false;

  // NOTE current not support
  easehts::VariantContext panel_of_normals_VC;
  easehts::VariantContext dbsnp_VC;

  bool covered = false;
  double power;
  double tumor_power;
  double normal_power;
  double normal_power_with_snp_prior;
  double normal_power_no_snp_prior;

  char alt_allele = 'N';
  std::string tumor_sample_name = "TUMOR";
  std::string normal_sample_name = "NORMAL";

  double contamination_fraction;

  double contaminant_lod;

  int tumor_q20_count;
  int normal_q20_count;

  int total_reads;
  int map_q0_reads;
  int initial_tumor_ref_counts;
  int initial_tumor_alt_counts;
  int initial_tumor_ref_quality_sum;
  int initial_tumor_alt_quality_sum;
  int initial_tumor_non_ref_quality_sum;
  int initial_tumor_read_depth;
  int initial_normal_ref_counts;
  int initial_normal_alt_counts;
  int initial_normal_ref_quality_sum;
  int initial_normal_alt_quality_sum;
  int tumor_ref_max_map_q;
  int tumor_alt_max_map_q;
  int initial_normal_read_depth;
  easehts::DiploidGenotype initial_normal_best_genotype;

  double initial_tumor_lod;
  double initial_normal_lod;

  double tumor_F;
  double tumor_lodF_star;
  double tumor_lodF_star_forward;
  double tumor_lodF_star_reverse;

  double normal_F;

  double power_to_detect_positive_strand_artifact;
  double power_to_detect_negative_strand_artifact;

  std::vector<int> strand_contingency_table;

  std::vector<int> tumor_alt_forward_offsets_in_read;
  std::vector<int> tumor_alt_reverse_offsets_in_read;

  double tumor_forward_offsets_in_read_median;
  double tumor_forward_offsets_in_read_mad;
  double normal_forward_offsets_in_read_median;
  double normal_forward_offsets_in_read_mad;

  int tumor_insertion_count;
  int tumor_deletion_count;

  std::vector<int> tumor_ref_quality_scores;
  std::vector<int> tumor_alt_quality_scores;
  std::vector<int> normal_ref_quality_scores;
  std::vector<int> normal_alt_quality_scores;

  std::vector<std::string> rejection_reasons;
  bool rejected = false;

};

} // mutect
} // ncic

#endif
