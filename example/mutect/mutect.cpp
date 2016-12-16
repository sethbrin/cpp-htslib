#include "mutect.h"
#include "locus_read_pile.h"
#include "candidate_mutation.h"
#include "variable_allelic_ratio_genotype_likelyhoods.h"
#include "cga_alignment_utils.h"
#include "sequence_utils.h"

#include <easehts/genome_loc.h>
#include <easehts/gatk/pileup.h>

#include <algorithm>
#include <atomic>
#include <array>
#include <cmath>
#include <functional>
#include <thread>
#include <vector>
#include <unordered_map>
#include <map>
#include <memory>

namespace ncic {
namespace mutect {

const std::string Worker::kValidBases = "ACGT";
const int Worker::kMinQSumQScore = 13;
const int Worker::kReferenceHalfWindowLength = 150;
const int Worker::kMaxReadMismatchQualityScoreSum = 100;
const char Worker::kMappedByMate = 'M';

void Worker::Run(const easehts::GenomeLoc& interval) {
  //printf("run: %s\n", interval.ToString().c_str());

  std::vector<easehts::gatk::GATKPileupTraverse> normal_traverses;
  std::vector<easehts::gatk::GATKPileupTraverse> tumor_traverses;

  bool is_downsampling = mutect_args_.downsampling.getValue();
  for (easehts::BAMIndexReader& reader : normal_readers_) {
    normal_traverses.emplace_back(&reader, interval,
                                  is_downsampling);
  }

  for (easehts::BAMIndexReader& reader : tumor_readers_) {
    tumor_traverses.emplace_back(&reader, interval,
                                 is_downsampling);
  }

  // the min postion of current site
  uint64_t cur_site = (uint64_t)-1;
  while (true) {
    uint64_t min_contig_pos = (uint64_t)-1;
    bool is_finished = true;
    for(auto& traverse : normal_traverses) {
      if (traverse.CurrentPileup().GetContigPos() == cur_site) {
        if (traverse.HasNext()) {
          is_finished = false;
        }
      }
      if (traverse.CurrentPileup().GetContigPos() < min_contig_pos) {
        min_contig_pos = traverse.CurrentPileup().GetContigPos();
      }
    }

    for(auto& traverse : tumor_traverses) {
      if (traverse.CurrentPileup().GetContigPos() == cur_site) {
        if (traverse.HasNext()) {
          is_finished = false;
        }
      }
      if (traverse.CurrentPileup().GetContigPos() < min_contig_pos) {
        min_contig_pos = traverse.CurrentPileup().GetContigPos();
      }
    }

    if (is_finished || min_contig_pos == (uint64_t)-1) {
      break;
    }
    cur_site = min_contig_pos;
    int cur_pos = (int)cur_site;
    int cur_contig_id = (int)(cur_site >> 32);

    ERROR_COND(cur_contig_id != interval.GetContigId(),
               easehts::utils::StringFormatCStr(
                   "the contig[%d] is not equal to the interval contig id[%d]",
                   cur_contig_id, interval.GetContigId()));

    if (cur_pos < interval.GetStart()) {
      continue;
    }
    if (cur_pos >= interval.GetStop()) {
      return;
    }

    easehts::GenomeLoc location(interval.GetContig(), interval.GetContigId(),
                                cur_pos, cur_pos);

    PrepareResult(location, min_contig_pos, tumor_traverses, normal_traverses);
  }
}

void Worker::PrepareVariantContext(const easehts::GenomeLoc& location) {
  for (auto& item : cosmic_traverses_) {
    item->SeekFroward(location);
    cosmic_variant_context_ = item->GetFirstRecord();
    if (cosmic_variant_context_) break;
  }

  for (auto& item : dbsnp_traverses_) {
    item->SeekFroward(location);
    dbsnp_variant_context_ = item->GetFirstRecord();
    if (dbsnp_variant_context_) break;
  }

}

void Worker::PrepareResult(
    const easehts::GenomeLoc& location,
    const uint64_t min_contig_pos,
    const std::vector<easehts::gatk::GATKPileupTraverse>& tumor_traverses,
    const std::vector<easehts::gatk::GATKPileupTraverse>& normal_traverses) {
  const char up_ref =
    std::toupper(reference_.GetSequenceAt(
            location.GetContig(),
            location.GetStart(), location.GetStop())[0]);
  // only process bases where the reference is [ACGT], because the FASTA
  // for HG18 has N,M and R!
  if (kValidBases.find(up_ref) == std::string::npos) {
    return;
  }
  LocusReadPile normal_read_pile(
      SampleType::NORMAL, up_ref, mutect_args_.min_qscore.getValue(),
      0, true, true, mutect_args_.enable_qscore_output.getValue());
  LocusReadPile tumor_read_pile(
      SampleType::TUMOR, up_ref, mutect_args_.min_qscore.getValue(),
      kMinQSumQScore, false, mutect_args_.artifact_detection_mode.getValue(),
      mutect_args_.enable_qscore_output.getValue());
  for(auto& traverse : tumor_traverses) {
    if (traverse.CurrentPileup().GetContigPos() == min_contig_pos) {
      tumor_read_pile.AddPileupElement(traverse.CurrentPileup());
    }
  }

  for(auto& traverse : normal_traverses) {
    if (traverse.CurrentPileup().GetContigPos() == min_contig_pos) {
      normal_read_pile.AddPileupElement(traverse.CurrentPileup());
    }
  }

  if (normal_read_pile.Size() == 0 &&
      tumor_read_pile.Size() == 0 &&
      !mutect_args_.force_output.getValue()) {
    return;
  }
  WARN_COND(false, easehts::utils::StringFormatCStr(
          "tid:%d pos:%d number of pileup:%d",
          location.GetContigId(), location.GetStart(),
          normal_read_pile.Size() +
          tumor_read_pile.Size()));

  tumor_read_pile.InitPileups();
  normal_read_pile.InitPileups();

  //printf("location: %d\n", location.GetStart() + 1);
  //for (int i=0; i < tumor_read_pile.pileup_.Size(); i++) {
  //  bam1_t* read = tumor_read_pile.pileup_[i].GetRead();
  //  printf("%s-%d\n", easehts::SAMBAMRecord::GetQueryName(read),
  //  easehts::SAMBAMRecord::GetSequenceLength(read));
  //}
  PrepareCondidate(up_ref, location, tumor_read_pile, normal_read_pile);
}

void Worker::PrepareCondidate(
    const char up_ref,
    const easehts::GenomeLoc& location,
    const LocusReadPile& tumor_read_pile,
    const LocusReadPile& normal_read_pile) {
  PrepareVariantContext(location);
  // remove the effect of cosmic from dbSNP
  // XXX current not support, just set false
  bool germline_at_risk = false;

  // compute coverage flags
  int tumor_covered_depth_threshold = 14;
  int normal_covered_depth_threshold = germline_at_risk ? 19 : 18;
  if (!has_normal_bam_) {
    normal_covered_depth_threshold = 0;
  }

  int tumor_base_count = tumor_read_pile.final_pileup_.Size();
  int normal_base_count = normal_read_pile.final_pileup_.Size();

  bool is_tumor_covered = tumor_base_count >= tumor_covered_depth_threshold;
  bool is_normal_covered = normal_base_count >= normal_covered_depth_threshold;
  bool is_base_covered = is_tumor_covered && is_normal_covered;
  if (!has_normal_bam_) {
    is_base_covered  = is_tumor_covered;
  }
  int tumor_q20_base_count = tumor_read_pile.final_pileup_
    .GetBaseFilteredPileupCount(20);
  int normal_q20_base_count = normal_read_pile.final_pileup_
    .GetBaseFilteredPileupCount(20);

  // calculate power
  double tumor_power = tumor_power_calculator_
    .CachingPowerCalculation(
        tumor_base_count,
        mutect_args_.power_contant_af.getValue());
  double normal_power_no_snp_prior = normal_novel_site_power_calculator_
    .CachingPowerCalculation(normal_base_count);
  double normal_power_with_snp_prior = normal_db_snp_site_power_calculator_
    .CachingPowerCalculation(normal_base_count);

  double normal_power = germline_at_risk ?
    normal_power_with_snp_prior : normal_power_no_snp_prior;

  double combine_power = tumor_power * normal_power;
  if (!has_normal_bam_) {
    combine_power = tumor_power;
  }

  int map_q0_reads = tumor_read_pile.quality_score_filter_pileup_
    .GetNumberofMappingQualityZeroReads();
  map_q0_reads += normal_read_pile.quality_score_filter_pileup_
    .GetNumberofMappingQualityZeroReads();


  int total_reads =
    tumor_read_pile.quality_score_filter_pileup_.Size() +
    normal_read_pile.quality_score_filter_pileup_.Size();

  bool force_output = mutect_args_.force_output.getValue();

  std::string sequence_context = reference_.CreateSequenceContext(location, 3);
  int window_start = location.GetStart() - Worker::kReferenceHalfWindowLength;
  int window_stop = location.GetStart() + Worker::kReferenceHalfWindowLength;
  easehts::GenomeLoc window(location.GetContig(), location.GetContigId(),
                            window_start, window_stop);
  easehts::ReferenceSequence ref_bases = reference_.GetSequenceAt(window);

  std::map<double, CandidateMutation*> message_by_tumor_lod;
  // Test each of the possible alternate alleles
  for (const char alt_allele : kValidBases) {
    if (alt_allele == up_ref) continue;
    if (!force_output &&
        tumor_read_pile.quality_sums_.GetCounts(alt_allele) == 0) continue;

    CandidateMutation* pCandidate = new CandidateMutation(location, up_ref);
    CandidateMutation& candidate = *pCandidate;
    candidate.sequence_context = sequence_context;
    candidate.covered = is_base_covered;
    candidate.power = combine_power;
    candidate.tumor_power = tumor_power;
    candidate.normal_power = normal_power;
    candidate.normal_power_with_snp_prior = normal_power_with_snp_prior;
    candidate.normal_power_no_snp_prior = normal_power_no_snp_prior;
    candidate.tumor_q20_count = tumor_q20_base_count;
    candidate.normal_q20_count = normal_q20_base_count;
    candidate.initial_tumor_non_ref_quality_sum =
      tumor_read_pile.quality_sums_.GetOtherQualities(up_ref);
    candidate.alt_allele = alt_allele;
    candidate.map_q0_reads = map_q0_reads;
    candidate.total_reads = total_reads;
    candidate.contamination_fraction =
      mutect_args_.fraction_contamination.getValue();
    // TODO
    // candidate.contamination_fraction
    // candidate.cosmic_site
    // candidate.dbsnp_site
    candidate.tumor_F = tumor_read_pile.EstimateAlleleFraction(
        up_ref, alt_allele);

    if (!force_output &&
        candidate.tumor_F < mutect_args_.tumor_f_pretest.getValue()) continue;

    candidate.initial_tumor_alt_counts =
      tumor_read_pile.quality_sums_.GetCounts(alt_allele);
    candidate.initial_tumor_ref_counts =
      tumor_read_pile.quality_sums_.GetCounts(up_ref);
    candidate.initial_tumor_alt_quality_sum =
      tumor_read_pile.quality_sums_.GetQualitySum(alt_allele);
    candidate.initial_tumor_ref_quality_sum =
      tumor_read_pile.quality_sums_.GetQualitySum(up_ref);

    double tumor_lod = tumor_read_pile.CalculateAltVsRefLOD(
        alt_allele, candidate.tumor_F, 0);
    candidate.tumor_lodF_star = tumor_lod;

    candidate.initial_tumor_read_depth = tumor_read_pile.final_pileup_.Size();
    candidate.tumor_insertion_count = tumor_read_pile.insertion_count_;
    candidate.tumor_deletion_count = tumor_read_pile.deletions_count_;

    if (candidate.tumor_lodF_star <
        mutect_args_.initial_tumor_lod_threshold.getValue()) continue;

    // calculate lod of contaminant
    double contaminant_F = std::min(contaminat_alternate_fraction_,
                                    candidate.tumor_F);
    VariableAllelicRatioGenotypeLikelihoods contaminant_likelihoods(
        up_ref, contaminant_F);

    std::vector<easehts::gatk::PileupElement*> pe_list =
      tumor_read_pile.final_pileup_.GetElements();
    auto pileup_comparator_by_alt_ref = [alt_allele](
        const easehts::gatk::PileupElement* lhs,
        const easehts::gatk::PileupElement* rhs)->bool const {
      if (lhs->GetBase() == rhs->GetBase()) {
        // if the bases are the same, the higher quality score comes first
        return lhs->GetQual() > rhs->GetQual();
      } else {
        if (lhs->GetBase() == alt_allele) {
          return true;
        } else if (rhs->GetBase() == alt_allele) {
          return false;
        } else {
          return lhs->GetBase() < rhs->GetBase();
        }
      }
    };
    std::sort(pe_list.begin(), pe_list.end(),
              pileup_comparator_by_alt_ref);

    int reads_to_keep = static_cast<int>(pe_list.size() *
                                         contaminat_alternate_fraction_);
    for (const auto& pe : pe_list) {
      char base = pe->GetBase();
      if (base == alt_allele) {
        // if we've retained all we need, then truen the
        // remainder of alts to ref
        if (reads_to_keep == 0) {
          base = up_ref;
        } else {
          reads_to_keep--;
        }
      }

      contaminant_likelihoods.Add(base, pe->GetQual());
    }
    std::array<double, 3> ref_het_hom =
      LocusReadPile::ExtractRefHemHom(contaminant_likelihoods,
                                      up_ref, alt_allele);
    candidate.contaminant_lod = ref_het_hom[1] - ref_het_hom[0];

    VariableAllelicRatioGenotypeLikelihoods normal_gl =
      normal_read_pile.CalculateLikelihoods(
          normal_read_pile.quality_score_filter_pileup_);
    candidate.initial_normal_best_genotype =
      normal_read_pile.GetBestGenotype(normal_gl);
    candidate.initial_normal_lod = LocusReadPile::GetRefVsAlt(
        normal_gl, up_ref, alt_allele);

    candidate.normal_F = std::max(LocusReadPile::EstimateAlleleFraction(
            normal_read_pile.quality_score_filter_pileup_, up_ref,
            alt_allele),
        (double)mutect_args_.minimum_normal_allele_fraction.getValue());


    const QualitySums& normal_qs = normal_read_pile.quality_sums_;
    candidate.initial_normal_alt_quality_sum =
      normal_qs.GetQualitySum(alt_allele);
    candidate.initial_normal_ref_quality_sum =
      normal_qs.GetQualitySum(up_ref);

    candidate.normal_alt_quality_scores =
      normal_qs.GetBaseQualityScores(alt_allele);
    candidate.normal_ref_quality_scores =
      normal_qs.GetBaseQualityScores(up_ref);

    candidate.initial_normal_alt_counts = normal_qs.GetCounts(alt_allele);
    candidate.initial_normal_ref_counts = normal_qs.GetCounts(up_ref);
    candidate.initial_normal_read_depth = normal_read_pile.final_pileup_.Size();


    // TODO: reduce the unnecessary compute pileup
    LocusReadPile t2(SampleType::NORMAL, up_ref, 0,
                     0, false, false,
                     mutect_args_.enable_qscore_output.getValue());
    FilterReads(ref_bases, location, window,
                tumor_read_pile.final_pileup_, true, &(t2.pileup_));
    t2.InitPileups();

    // if there are no reads remaining, abandon this theory
    if (!force_output && t2.final_pileup_.Size() == 0) continue;

    candidate.initial_tumor_alt_counts = t2.quality_sums_.GetCounts(alt_allele);
    candidate.initial_tumor_ref_counts = t2.quality_sums_.GetCounts(up_ref);
    candidate.initial_tumor_alt_quality_sum =
      t2.quality_sums_.GetQualitySum(alt_allele);
    candidate.initial_tumor_ref_quality_sum =
      t2.quality_sums_.GetQualitySum(up_ref);

    candidate.tumor_alt_quality_scores =
      t2.quality_sums_.GetBaseQualityScores(alt_allele);
    candidate.tumor_ref_quality_scores =
      t2.quality_sums_.GetBaseQualityScores(up_ref);

    VariableAllelicRatioGenotypeLikelihoods t2_gl =
      t2.CalculateLikelihoods(t2.final_pileup_);
    candidate.initial_tumor_lod = t2.GetAltVsRef(t2_gl, up_ref, alt_allele);
    candidate.initial_tumor_read_depth = t2.final_pileup_.Size();

    candidate.tumor_F = t2.EstimateAlleleFraction(up_ref, alt_allele);
    candidate.tumor_lodF_star = t2.CalculateAltVsRefLOD(alt_allele,
                                                        candidate.tumor_F, 0);

    // TODO: clean up use of forward/reverse vs positive/negative
    // (prefer the latter)
    easehts::gatk::ReadBackedPileup tmp_pileup;

    easehts::gatk::ReadBackedPileup forward_pileup;
    FilterReads(ref_bases, location, window,
                tumor_read_pile.final_pileup_positive_strand_,
                true, &(tmp_pileup));
    tmp_pileup.GetPileupByAndFilter(
        &forward_pileup,
        {easehts::gatk::PileupFilter::IsNotDeletion,
        std::bind(easehts::gatk::PileupFilter::IsBaseQualityLarge,
                  std::placeholders::_1, 0),
        easehts::gatk::PileupFilter::IsMappingQualityLargerThanZero,
        easehts::gatk::PileupFilter::IsPositiveStrand});
    double f2_forward = LocusReadPile::EstimateAlleleFraction(
        forward_pileup, up_ref, alt_allele);
    candidate.tumor_lodF_star_forward = t2.CalculateAltVsRefLOD(
        forward_pileup, alt_allele, f2_forward, 0.0);

    tmp_pileup.Clear();
    easehts::gatk::ReadBackedPileup reverse_pileup;
    FilterReads(ref_bases, location, window,
                tumor_read_pile.final_pileup_negative_strand_,
                true, &(tmp_pileup));
    tmp_pileup.GetPileupByAndFilter(
        &reverse_pileup,
        {easehts::gatk::PileupFilter::IsNotDeletion,
        std::bind(easehts::gatk::PileupFilter::IsBaseQualityLarge,
                  std::placeholders::_1, 0),
        easehts::gatk::PileupFilter::IsMappingQualityLargerThanZero,
        easehts::gatk::PileupFilter::IsNegativeStrand});
    double f2_reverse = LocusReadPile::EstimateAlleleFraction(
        reverse_pileup, up_ref, alt_allele);
    candidate.tumor_lodF_star_reverse = t2.CalculateAltVsRefLOD(
        reverse_pileup, alt_allele, f2_reverse, 0.0);


    candidate.power_to_detect_positive_strand_artifact =
      strand_artifact_power_calculator_.CachingPowerCalculation(
          reverse_pileup.Size(), candidate.tumor_F);
    candidate.power_to_detect_negative_strand_artifact =
      strand_artifact_power_calculator_.CachingPowerCalculation(
          forward_pileup.Size(), candidate.tumor_F);

    candidate.strand_contingency_table =
      SequenceUtils::GetStrandContingencyTable(forward_pileup, reverse_pileup,
                                               up_ref, alt_allele);


    easehts::gatk::ReadBackedPileup mutant_pileup;
    easehts::gatk::ReadBackedPileup reference_pileup;
    for (int idx=0; idx<t2.final_pileup_.Size(); idx++) {
      easehts::gatk::PileupElement* p = t2.final_pileup_[idx];
      easehts::SAMBAMRecord* read = p->GetRead();
      int offset = p->GetOffset();
      char cur_char = read->GetSequenceAt(offset);

      if (cur_char == alt_allele) {
        mutant_pileup.AddElement(p);
      } else if (cur_char == up_ref) {
        reference_pileup.AddElement(p);
      }

    }
    // start with just the tumor pile
    candidate.tumor_ref_max_map_q = reference_pileup.GetMaxMappingQuals();
    candidate.tumor_alt_max_map_q = mutant_pileup.GetMaxMappingQuals();

    // Set the maximum observed mapping quality score for the reference and
    // alternate alleles
    candidate.tumor_alt_forward_offsets_in_read =
      SequenceUtils::GetForwardOffsetsInRead(mutant_pileup, location);
    candidate.tumor_alt_reverse_offsets_in_read =
      SequenceUtils::GetReverseOffsetsInRead(mutant_pileup, location);

    if (candidate.tumor_alt_forward_offsets_in_read.size() > 0) {
      std::vector<int> offsets = candidate.tumor_alt_forward_offsets_in_read;
      double median = MutectStats::GetMedian<int>(offsets);
      candidate.tumor_forward_offsets_in_read_median = median;
      candidate.tumor_forward_offsets_in_read_mad =
        MutectStats::CalculateMAD(offsets, median);
    }

    if (candidate.tumor_alt_reverse_offsets_in_read.size() > 0) {
      std::vector<int> offsets = candidate.tumor_alt_reverse_offsets_in_read;
      double median = MutectStats::GetMedian<int>(offsets);
      candidate.tumor_reverse_offsets_in_read_median = median;
      candidate.tumor_reverse_offsets_in_read_mad =
        MutectStats::CalculateMAD(offsets, median);
    }

    // test to see if the candidate should be rejected
    PerformRejection(candidate);

    if (mutect_args_.force_alleles.getValue()) {
      call_stats_generator_.WriteCallStats(candidate);
    } else {
      message_by_tumor_lod[candidate.initial_tumor_lod] = &candidate;
    }
  }

  // if more than one site passes the tumor lod threshold for KEEP the fail the
  // tri_allelic site filter
  int passing_candidates = 0;
  for (const auto& c : message_by_tumor_lod) {
    if (c.second->tumor_lodF_star >=
        mutect_args_.tumor_lod_threshold.getValue()) {
      passing_candidates ++;
    }
  }

  if (passing_candidates > 1) {
    for (const auto& c : message_by_tumor_lod) {
      c.second->AddRejectionReason("triallelic_site");
    }
  }

  // write out call stats for the "best" candidate
  if (!message_by_tumor_lod.empty()) {
    const CandidateMutation& m = *(message_by_tumor_lod.rbegin()->second);

    // only output passing calls OR rejected sites if ONLY_PASSING_CALLS is not
    // specified
    if (!m.IsRejected() ||
        (!mutect_args_.only_passing_calls.getValue())) {
      call_stats_generator_.WriteCallStats(m);
    }
  }

  for (const auto& c : message_by_tumor_lod) {
    delete c.second;
  }
}

void Worker::PerformRejection(
    CandidateMutation& candidate) {
  if (candidate.tumor_lodF_star <
      mutect_args_.tumor_lod_threshold.getValue()) {
    candidate.AddRejectionReason("fstar_tumor_lod");
  }

  if (mutect_args_.artifact_detection_mode.getValue()) return;

  if (candidate.tumor_insertion_count >=
      mutect_args_.gap_events_threshold.getValue() ||
      candidate.tumor_deletion_count >=
      mutect_args_.gap_events_threshold.getValue()) {
    candidate.AddRejectionReason("nearby_gap_events");
  }

  if (mutect_args_.fraction_contamination.getValue() +
      mutect_args_.minimum_mutation_cell_fraction.getValue() > 0 &&
      candidate.tumor_lodF_star <= mutect_args_.tumor_lod_threshold.getValue() +
      std::max(0.0, candidate.contaminant_lod)) {
    candidate.AddRejectionReason("possible_contamination");
  }

  if (candidate.IsGermlineAtRisk() &&
      candidate.initial_normal_lod <
      mutect_args_.normal_dbsnp_lod_threshold.getValue()) {
    candidate.AddRejectionReason("germline_risk");
  }

  if (candidate.initial_normal_lod <
      mutect_args_.normal_lod_threshold.getValue()) {
    candidate.AddRejectionReason("normal_lod");
  }

  if ((candidate.initial_normal_alt_counts >=
       mutect_args_.max_alt_alleles_in_normal_count.getValue() ||
      candidate.normal_F >=
      mutect_args_.max_alt_allele_in_normal_fraction.getValue()) &&
      candidate.initial_normal_alt_quality_sum >
      mutect_args_.max_alt_alleles_in_normal_qscore_sum.getValue()) {
    candidate.AddRejectionReason("alt_allele_in_normal");
  }

  if ((candidate.tumor_forward_offsets_in_read_median <=
       mutect_args_.pir_median_threshold.getValue() &&
       candidate.tumor_forward_offsets_in_read_mad <=
       mutect_args_.pir_mad_threshold.getValue()) ||
      (candidate.tumor_reverse_offsets_in_read_median <=
       mutect_args_.pir_median_threshold.getValue() &&
       candidate.tumor_reverse_offsets_in_read_mad <=
       mutect_args_.pir_mad_threshold.getValue())) {
    candidate.AddRejectionReason("clustered_read_position");
  }

  // TODO: sync naming(is it positive or forward)
  if ((candidate.power_to_detect_negative_strand_artifact >=
       mutect_args_.strand_artifact_power_threshold.getValue() &&
      candidate.tumor_lodF_star_forward <
      mutect_args_.strand_artifact_lod_threshold.getValue()) ||
      (candidate.power_to_detect_positive_strand_artifact >=
       mutect_args_.strand_artifact_power_threshold.getValue() &&
      candidate.tumor_lodF_star_reverse <
      mutect_args_.strand_artifact_lod_threshold.getValue())) {
    candidate.AddRejectionReason("strand_artifact");
  }

  if (candidate.total_reads > 0 &&
      static_cast<float>(candidate.map_q0_reads)/candidate.total_reads >=
      mutect_args_.fraction_mapq0_threshold.getValue()) {
    candidate.AddRejectionReason("poor_mapping_region_mapq0");
  }

  if (candidate.tumor_alt_max_map_q <
      mutect_args_.required_maximum_alt_allele_mapping_quality_score.getValue()) {
    candidate.AddRejectionReason("poor_mapping_region_alternate_allele_mapq");
  }

  if (candidate.IsSeenInPanelOfNormals()) {
    if (!candidate.cosmic_site) {
      candidate.AddRejectionReason("seen_in_panel_of_normals");
    }
  }

}

int Worker::Input(void *data, bam1_t *b) {
  InputData& input_data = *(InputData*)data;
  while (true) {
    if (!input_data.reader->HasNext(b)) {
      return -1;
    }
    // check if b is valid, recording to Mutect project LocusWalker class
    easehts::SAMBAMRecord record(b);
    if (record.GetReadUnmappedFlag() ||
        record.GetAlignmentStart() ==
        easehts::SAMBAMRecord::NO_ALIGNMENT_START ||
        record.GetNotPrimaryAlignmentFlag() ||
        record.GetDuplicateReadFlag() ||
        record.GetReadFailsVendorQualityCheckFlag()) {
      // read next
    } else {
      // FIXME bug design, should set record to null, otherwise will delete the
      // memory
      record.SetRawRecord(nullptr);
      return 0;
    }
  }
  return -1;
}

void Worker::FilterReads(
    const easehts::ReferenceSequence& ref_bases,
    const easehts::GenomeLoc& window,
    const easehts::GenomeLoc& location,
    const easehts::gatk::ReadBackedPileup pileup,
    bool filter_mate_rescue_reads,
    easehts::gatk::ReadBackedPileup* pPileup) {
  int size = pileup.Size();
  for (int idx = 0; idx < size; idx++) {
    easehts::gatk::PileupElement* p = pileup[idx];
    easehts::SAMBAMRecord* read = p->GetRead();

    int mismatch_quality_sum = CGAAlignmentUtils::MismatchesInRefWindow(
        ref_bases, p, window, location, false, true);

    // do we have too many mismatch overall?
    if (mismatch_quality_sum > kMaxReadMismatchQualityScoreSum) continue;

    // is this a heavily cliped read?
    if (SequenceUtils::IsReadHeavilySoftClipped(
            read, mutect_args_.heavily_clipped_read_fraction.getValue())) {
      continue;
    }

    // was this read only placed because it's mate was uniquely placed?
    // (supplied by BWA)
    // TODO: current not support SAMTags
    pPileup->AddElement(p);
  }
}

void Mutect::Run() {
  easehts::GenomeLocParser parser(reference_.GetSequenceDictionary());
  std::vector<easehts::GenomeLoc> intervals =
    easehts::IntervalUtils::LoadIntervals(
        parser, mutect_args_.interval_file.getValue(),
        mutect_args_.interval_padding.getValue());

  call_stats_generator_.WriteHeader();
  int thread_cnt = mutect_args_.thread_cnt.getValue();

  std::atomic<int> interval_index(0);
  std::vector<std::thread> workers;
  workers.reserve(thread_cnt);
  for (int i = 0; i < thread_cnt; i++) {
    workers.emplace_back([this, &intervals, &interval_index]() {
      Worker worker(this->mutect_args_, this->reference_,
                    this->call_stats_generator_);
      while (interval_index < intervals.size()) {
        size_t index = interval_index++;
        fprintf(stderr, "%d--%d-%d\n", index, intervals[index].GetStart(),
                intervals[index].GetStop());
        if (index >= intervals.size()) break;

        // As the htslib file reader is [start, end)
        // but the interval is [start, end]
        worker.Run(easehts::GenomeLoc(intervals[index].GetContig(),
                                      intervals[index].GetContigId(),
                                      intervals[index].GetStart(),
                                      intervals[index].GetStop() + 1));
      }
    });
  }

  for (int i = 0; i < thread_cnt; i++) {
    workers[i].join();
  }
}


} // mutect
} // ncic
