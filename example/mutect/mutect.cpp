#include "mutect.h"
#include "locus_read_pile.h"
#include "candidate_mutation.h"

#include <easehts/genome_loc.h>

#include <atomic>
#include <cmath>
#include <thread>
#include <vector>

namespace ncic {
namespace mutect {

const std::string Worker::kValidBases = "ACGT";
const int Worker::kMinQSumQScore = 13;

void Worker::Run(const easehts::GenomeLoc& interval) {
  //printf("run: %s\n", interval.ToString().c_str());

  std::vector<easehts::PileupTraverse> normal_traverses;
  std::vector<easehts::PileupTraverse> tumor_traverses;

  int overlap_size = mutect_args_.overlap_size.getValue();
  for (easehts::BAMIndexReader& reader : normal_readers_) {
    InputData input_data;
    reader.SetRegion(interval.GetContigId(),
                     interval.GetStart() - overlap_size,
                     interval.GetStop() + overlap_size);
    input_data.reader = &reader;
    normal_traverses.emplace_back(&Worker::Input, &input_data);
  }

  for (easehts::BAMIndexReader& reader : tumor_readers_) {
    InputData input_data;
    reader.SetRegion(interval.GetContigId(),
                     interval.GetStart() - overlap_size,
                     interval.GetStop() + overlap_size);
    input_data.reader = &reader;
    tumor_traverses.emplace_back(&Worker::Input, &input_data);
  }

  // the min postion of current site
  uint64_t cur_site = (uint64_t)-1;
  while (true) {
    uint64_t min_contig_pos = (uint64_t)-1;
    for(auto& traverse : normal_traverses) {
      if (traverse.CurrentPileup().GetContigPos() == cur_site) {
        traverse.GetNext();
      }
      if (traverse.CurrentPileup().GetContigPos() < min_contig_pos) {
        min_contig_pos = traverse.CurrentPileup().GetContigPos();
      }
    }

    for(auto& traverse : tumor_traverses) {
      if (traverse.CurrentPileup().GetContigPos() == cur_site) {
        traverse.GetNext();
      }
      if (traverse.CurrentPileup().GetContigPos() < min_contig_pos) {
        min_contig_pos = traverse.CurrentPileup().GetContigPos();
      }
    }

    if (min_contig_pos == (uint64_t)-1) {
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
    if (cur_pos > interval.GetStop()) {
      return;
    }

    easehts::GenomeLoc location(interval.GetContig(), interval.GetContigId(),
                                cur_pos, cur_pos);

    PrepareResult(location, min_contig_pos, tumor_traverses, normal_traverses);
  }
}

void Worker::PrepareResult(const easehts::GenomeLoc& location,
                           const uint64_t min_contig_pos,
                           const std::vector<easehts::PileupTraverse>& tumor_traverses,
                           const std::vector<easehts::PileupTraverse>& normal_traverses) {
  const char up_ref = std::toupper(reference_.GetSequenceAt(location.GetContig(),
      location.GetStart(), location.GetStop())[0]);
  // only process bases where the reference is [ACGT], because the FASTA
  // for HG18 has N,M and R!
  if (kValidBases.find(up_ref) == std::string::npos) {
    return;
  }
  LocusReadPile normal_read_pile(SampleType::NORMAL, up_ref, mutect_args_.min_qscore.getValue(),
                                 0, true, true, mutect_args_.enable_qscore_output.getValue());
  LocusReadPile tumor_read_pile(SampleType::TUMOR, up_ref, mutect_args_.min_qscore.getValue(),
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
  WARN(easehts::utils::StringFormatCStr("tid:%d pos:%d number of pileup:%d",
                                        location.GetContigId(), location.GetStart(),
                                        normal_read_pile.Size() +
                                        tumor_read_pile.Size()));

  tumor_read_pile.InitPileups();
  normal_read_pile.InitPileups();

  PrepareCondidate(up_ref, location, tumor_read_pile, normal_read_pile);
}

void Worker::PrepareCondidate(
    const char up_ref,
    const easehts::GenomeLoc& location,
    const LocusReadPile& tumor_read_pile,
    const LocusReadPile& normal_read_pile) {
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

  double normal_power = germline_at_risk ? normal_power_with_snp_prior : normal_power_no_snp_prior;

  double combine_power = tumor_power * normal_power;
  if (!has_normal_bam_) {
    combine_power = tumor_power;
  }

  int map_q0_reads =
    tumor_read_pile.quality_score_filter_pileup_.GetNumberofMappingQualityZeroReads() +
    normal_read_pile.quality_score_filter_pileup_.GetNumberofMappingQualityZeroReads();

  int total_reads =
    tumor_read_pile.quality_score_filter_pileup_.Size() +
    normal_read_pile.quality_score_filter_pileup_.Size();

  std::string sequence_context = reference_.CreateSequenceContext(location, 3);
  // Test each of the possible alternate alleles
  for (const char alt_allele : kValidBases) {
    if (alt_allele == up_ref) continue;
    if (!mutect_args_.force_output.getValue() &&
        tumor_read_pile.quality_sums_.GetCounts(alt_allele) == 0) continue;

    CandidateMutation candidate(location, up_ref);
    candidate.sequence_context = sequence_context;
    candidate.covered = is_base_covered;
    candidate.power = combine_power;
    candidate.tumor_power = tumor_power;
    candidate.normal_power = normal_power;
    candidate.normal_power_with_snp_prior = normal_power_with_snp_prior;
    candidate.normal_power_no_snp_prior = normal_power_no_snp_prior;
    candidate.tumor_q20_count = tumor_q20_base_count;
    candidate.normal_q20_count = normal_q20_base_count;
    candidate.initial_tumor_non_ref_quality_sum = tumor_read_pile.quality_sums_.GetOtherQualities(up_ref);
    candidate.alt_allele = alt_allele;
    candidate.map_q0_reads = map_q0_reads;
    candidate.total_reads = total_reads;
    candidate.contamination_fraction = mutect_args_.fraction_contamination.getValue();
    // TODO
    // candidate.contamination_fraction
    // candidate.cosmic_site
    // candidate.dbsnp_site
    candidate.tumor_F = tumor_read_pile.EstimateAlleleFraction(up_ref, alt_allele);

    if (!mutect_args_.force_output.getValue() &&
        candidate.tumor_F < mutect_args_.tumor_f_pretest.getValue()) continue;

    candidate.initial_tumor_alt_counts = tumor_read_pile.quality_sums_.GetCounts(alt_allele);
    candidate.initial_tumor_ref_counts = tumor_read_pile.quality_sums_.GetCounts(up_ref);
    candidate.initial_tumor_alt_quality_sum = tumor_read_pile.quality_sums_.GetQualitySum(alt_allele);
    candidate.initial_tumor_ref_quality_sum = tumor_read_pile.quality_sums_.GetQualitySum(up_ref);

    double tumor_lod = tumor_read_pile.CalculateAltVsRefLOD(alt_allele, candidate.tumor_F, 0);
    candidate.tumor_lodF_star = tumor_lod;

    candidate.initial_tumor_read_depth = tumor_read_pile.final_pileup_.Size();
    candidate.tumor_insertion_count = tumor_read_pile.insertion_count_;
    candidate.tumor_deletion_count = tumor_read_pile.deletions_count_;

    if (candidate.tumor_lodF_star < mutect_args_.initial_tumor_lod_threshold.getValue()) continue;

    // calculate lod of contaminant
    double contaminant_F = std::min(contaminat_alternate_fraction_,
                                    candidate.tumor_F);
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
        record.GetAlignmentStart() == easehts::SAMBAMRecord::NO_ALIGNMENT_START ||
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

void Mutect::Run() {
  easehts::GenomeLocParser parser(reference_.GetSequenceDictionary());
  std::vector<easehts::GenomeLoc> intervals =
    easehts::IntervalUtils::IntervalFileToList(parser,
                                               mutect_args_.interval_file.getValue());

  int thread_cnt = mutect_args_.thread_cnt.getValue();

  std::atomic<int> interval_index(0);
  std::vector<std::thread> workers;
  workers.reserve(thread_cnt);
  for (int i = 0; i < thread_cnt; i++) {
    workers.emplace_back([this, &intervals, &interval_index]() {
      Worker worker(this->mutect_args_, this->reference_);
      while (interval_index < intervals.size()) {
        size_t index = interval_index++;
        if (index >= intervals.size()) break;

        worker.Run(intervals[index]);
      }
    });
  }

  for (int i = 0; i < thread_cnt; i++) {
    workers[i].join();
  }
}


} // mutect
} // ncic
