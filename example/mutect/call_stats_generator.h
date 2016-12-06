//
// Created by zp on 11/27/16.
//

#ifndef MUTECT_CALL_STATS_GENERATOR_H_
#define MUTECT_CALL_STATS_GENERATOR_H_

#include "candidate_mutation.h"

#include <easehts/noncopyable.h>
#include <easehts/utils.h>

#include <string>
#include <unordered_map>
#include <vector>
#include <stdio.h>

namespace ncic {
namespace mutect {

class CallStatsGenerator : public easehts::NonCopyable {
 public:
  CallStatsGenerator(bool enable_qscore_output, const std::string& file_name) {
    SetupOutputColumns(enable_qscore_output);

    fp_ = fopen(file_name.c_str(), "w");
    if (fp_ == nullptr) {
      ERROR(easehts::utils::StringFormatCStr("Can not open filename: %s", file_name.c_str()));
    }
  }

  ~CallStatsGenerator() {
    fclose(fp_);
  }

  void WriteHeader() {
    fprintf(fp_, "%s\n", GenerateHeader().c_str());
  }

  void WriteCallStats(const CandidateMutation& candidate) {
    fprintf(fp_, "%s\n", GenerateCallStats(candidate).c_str());
  }

 private:
  std::string GenerateHeader() {
    return easehts::utils::Join(headers_, kTab);
  }

  std::string GenerateCallStats(const CandidateMutation& candidate) {
    std::unordered_map<std::string, std::string> d;

    std::string keep_string = "REJECT";
    if (!candidate.IsRejected()) {
      keep_string = "KEEP";
    }

    std::string site_info = GetSiteInfoString(candidate.dbsnp_site, candidate.cosmic_site);

    std::string strand_info = GetStrandTableString(candidate.strand_contingency_table);

    d["contig"] = candidate.location.GetContig();
    // XXX change to 1-based
    d["position"] = Format(candidate.location.GetStart() + 1);
    d["context"] = candidate.sequence_context;
    d["ref_allele"] = Format(candidate.ref_allele);
    d["alt_allele"] = Format(candidate.alt_allele);
    d["tumor_name"] = candidate.tumor_sample_name;
    d["normal_name"] = candidate.normal_sample_name;
    d["score"] = Format(candidate.score);
    d["dbsnp_site"] = site_info;
    d["covered"] = (candidate.covered?"COVERED":"UNCOVERED");
    d["power"] = Format(candidate.power);
    d["tumor_power"] = Format(candidate.tumor_power);
    d["normal_power"] = Format(candidate.normal_power);
    d["normal_power_nsp"] = Format(candidate.normal_power_no_snp_prior);
    d["normal_power_wsp"] = Format(candidate.normal_power_with_snp_prior);
    d["total_reads"] = Format(candidate.total_reads);
    d["map_Q0_reads"] = Format(candidate.map_q0_reads);
    d["init_t_lod"] = Format(candidate.initial_tumor_lod);
    d["t_lod_fstar"] = Format(candidate.tumor_lodF_star);
    d["t_lod_fstar_forward"] = Format(candidate.tumor_lodF_star_forward);
    d["t_lod_fstar_reverse"] = Format(candidate.tumor_lodF_star_reverse);
    d["tumor_f"] = Format(candidate.tumor_F);
    d["contaminant_fraction"] = Format(candidate.contamination_fraction);
    d["contaminant_lod"] = Format(candidate.contaminant_lod);
    d["t_q20_count"] = Format(candidate.tumor_q20_count);
    d["t_ref_count"] = Format(candidate.initial_tumor_ref_counts);
    d["t_alt_count"] = Format(candidate.initial_tumor_alt_counts);
    d["t_ref_sum"] = Format(candidate.initial_tumor_ref_quality_sum);
    d["t_alt_sum"] = Format(candidate.initial_tumor_alt_quality_sum);
    d["t_ref_max_mapq"] = Format(candidate.tumor_ref_max_map_q);
    d["t_alt_max_mapq"] = Format(candidate.tumor_alt_max_map_q);
    d["t_ins_count"] = Format(candidate.tumor_insertion_count);
    d["t_del_count"] = Format(candidate.tumor_deletion_count);
    d["normal_best_gt"] = Format(candidate.initial_normal_best_genotype.ToString());
    d["init_n_lod"] = Format(candidate.initial_normal_lod);
    d["normal_f"] = Format(candidate.normal_F);
    d["n_q20_count"] = Format(candidate.normal_q20_count);
    d["n_ref_count"] = Format(candidate.initial_normal_ref_counts);
    d["n_alt_count"] = Format(candidate.initial_normal_alt_counts);
    d["n_ref_sum"] = Format(candidate.initial_normal_ref_quality_sum);
    d["n_alt_sum"] = Format(candidate.initial_normal_alt_quality_sum);
    d["power_to_detect_positive_strand_artifact"] = Format(candidate.power_to_detect_positive_strand_artifact);
    d["power_to_detect_negative_strand_artifact"] = Format(candidate.power_to_detect_negative_strand_artifact);
    d["strand_bias_counts"] = Format(strand_info);
    d["tumor_alt_fpir_median"] = Format(candidate.tumor_forward_offsets_in_read_median);
    d["tumor_alt_fpir_mad"] = Format(candidate.tumor_forward_offsets_in_read_mad);
    d["tumor_alt_rpir_median"] = Format(candidate.tumor_reverse_offsets_in_read_median);
    d["tumor_alt_rpir_mad"] = Format(candidate.tumor_reverse_offsets_in_read_mad);
    d["observed_in_normals_count"] = Format(candidate.GetCountOfNormalsObservedIn());
    d["tumor_ref_qscores"] = Format(candidate.tumor_ref_quality_scores);
    d["tumor_alt_qscores"] = Format(candidate.tumor_alt_quality_scores);
    d["normal_ref_qscores"] = Format(candidate.normal_ref_quality_scores);
    d["normal_alt_qscores"] = Format(candidate.normal_alt_quality_scores);
    d["failure_reasons"] = easehts::utils::Join(candidate.rejection_reasons, kTab);
    d["judgement"] = keep_string;

    return Generate(d);
  }

  static const std::string kTab;

  void SetupOutputColumns(bool enable_qscore_output) {
    headers_.push_back("contig");
    headers_.push_back("position");
    headers_.push_back("context");
    headers_.push_back("ref_allele");
    headers_.push_back("alt_allele");
    headers_.push_back("tumor_name");
    headers_.push_back("normal_name");
    headers_.push_back("score");
    headers_.push_back("dbsnp_site");
    headers_.push_back("covered");
    headers_.push_back("power");
    headers_.push_back("tumor_power");
    headers_.push_back("normal_power");
    headers_.push_back("normal_power_nsp");
    headers_.push_back("normal_power_wsp");
    headers_.push_back("total_reads");
    headers_.push_back("map_Q0_reads");
    headers_.push_back("init_t_lod");
    headers_.push_back("t_lod_fstar");
    headers_.push_back("t_lod_fstar_forward");
    headers_.push_back("t_lod_fstar_reverse");
    headers_.push_back("tumor_f");
    headers_.push_back("contaminant_fraction");
    headers_.push_back("contaminant_lod");
    headers_.push_back("t_q20_count");
    headers_.push_back("t_ref_count");
    headers_.push_back("t_alt_count");
    headers_.push_back("t_ref_sum");
    headers_.push_back("t_alt_sum");
    headers_.push_back("t_ref_max_mapq");
    headers_.push_back("t_alt_max_mapq");
    headers_.push_back("t_ins_count");
    headers_.push_back("t_del_count");
    headers_.push_back("normal_best_gt");
    headers_.push_back("init_n_lod");
    headers_.push_back("normal_f");
    headers_.push_back("n_q20_count");
    headers_.push_back("n_ref_count");
    headers_.push_back("n_alt_count");
    headers_.push_back("n_ref_sum");
    headers_.push_back("n_alt_sum");
    headers_.push_back("power_to_detect_positive_strand_artifact");
    headers_.push_back("power_to_detect_negative_strand_artifact");
    headers_.push_back("strand_bias_counts");
    headers_.push_back("tumor_alt_fpir_median");
    headers_.push_back("tumor_alt_fpir_mad");
    headers_.push_back("tumor_alt_rpir_median");
    headers_.push_back("tumor_alt_rpir_mad");
    headers_.push_back("observed_in_normals_count");

    if (enable_qscore_output) {
      headers_.push_back("tumor_ref_qscores");
      headers_.push_back("tumor_alt_qscores");
      headers_.push_back("normal_ref_qscores");
      headers_.push_back("normal_alt_qscores");
    }

    headers_.push_back("failure_reasons");
    headers_.push_back("judgement");
  }

  std::string Generate(const std::unordered_map<std::string, std::string>& d) {
    std::vector<std::string> msg;
    msg.reserve(headers_.size());

    for (const auto& value : headers_) {
      msg.push_back(d.at(value));
    }
    fprintf(stderr, "*****position:%s\n", d.at("position").c_str());

    return easehts::utils::Join(msg, kTab);
  }

  std::string GetSiteInfoString(bool is_dbsnp_site, bool is_comic_site) {
    std::string site_info = "NOVEL";
    if (is_dbsnp_site) {
      site_info = "DBSNP";
    }
    if (is_comic_site) {
      site_info = "COSMIC";
    }

    if (is_dbsnp_site && is_comic_site) {
      site_info = "DBSNP+COSMIC";
    }
    return site_info;
  }

  std::string GetStrandTableString(std::array<int, 4> ci) {
    std::string res = "(";
    res += (std::to_string(ci[0]) + ",");
    res += (std::to_string(ci[1]) + ",");
    res += (std::to_string(ci[2]) + ",");
    res += (std::to_string(ci[3]) + ")");
    return res;
  }

  std::string Format(const std::string& s) {
    return s;
  }

  std::string Format(int i) {
    return std::to_string(i);
  }

  std::string Format(char ch) {
    return std::string(1, ch);
  }

  std::string Format(const std::vector<int>& ints) {
    if (ints.empty()) {
      return "n/a";
    }

    return easehts::utils::Join(ints, kTab);
  }

  std::string Format(double d) {
    if (d == CandidateMutation::kDoubleUnintialized) {
      return "n/a";
    }
    return std::to_string(d);
  }

  std::vector<std::string> headers_;
  FILE* fp_;
};


} // mutect
} // ncic

#endif
