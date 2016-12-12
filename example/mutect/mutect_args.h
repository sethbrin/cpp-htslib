//
// Created by zp on 11/20/16.
//

#ifndef MUTECT_MUTECE_ARGS_H_
#define MUTECT_MUTECE_ARGS_H_

#include <easehts/args_parser.h>
#include <easehts/noncopyable.h>

namespace ncic {
namespace mutect {

class MutectArgs : public ncic::easehts::NonCopyable {
 public:
  MutectArgs(int argc, char** argv)
    : interval_file('L', "intervals",   true , "interval file"),
    enable_qscore_output('e', "enable_qscore_output", false, "output quality scores of all bases used for ref and alt"),
    artifact_detection_mode('a', "artifact_detection_mode", false, "used when running the caller on a normal (as if it were a tumor) to detect artifacts"),
    force_output('f', "force_output", false, "force output for each site"),
    force_alleles('F', "force_alleles", false, "force output for all alleles at each site"),
    only_passing_calls(easehts::BaseOption::NO_OPTION, "only_passing_calls", false, "only emit passing calls"),
    initial_tumor_lod_threshold(easehts::BaseOption::NO_OPTION, "initial_tumor_lod", false, 4.0f, "Initial LOD threshold for calling tumor variant"),
    tumor_lod_threshold(easehts::BaseOption::NO_OPTION, "tumor_lod", false, 6.3f, "LOD threshold for calling tumor variant"),
    fraction_contamination(easehts::BaseOption::NO_OPTION, "fraction_contamination", false, 0.02f, "estimate of fraction (0-1) of physical contamination with other unrelated samples"),
    minimum_mutation_cell_fraction(easehts::BaseOption::NO_OPTION, "minimum_mutation_cell_fraction", false, 0.00f, "minimum fraction of cells which are presumed to have a mutation, used to handle non-clonality and contamination"),
    normal_lod_threshold(easehts::BaseOption::NO_OPTION, "normal_lod_threshold", false, 2.2f, "LOD threshold for calling normal non-germline"),
    normal_artifact_lod_threshold(easehts::BaseOption::NO_OPTION, "normal_artifact_lod", false, 1.0f, "LOD threshold for calling normal non-variant"),
    strand_artifact_lod_threshold(easehts::BaseOption::NO_OPTION, "strand_artifact_lod", false, 2.0f, "LOD threshold for calling strand bias"),
    strand_artifact_power_threshold(easehts::BaseOption::NO_OPTION, "strand_artifact_power", false, 0.9f, "power threshold for calling strand bias"),
    normal_dbsnp_lod_threshold(easehts::BaseOption::NO_OPTION, "normal_dbsnp_lod", false, 5.5f, "LOD threshold for calling normal non-variant at dbsnp sites"),
    minimum_normal_allele_fraction(easehts::BaseOption::NO_OPTION, "minimum_normal_allele_fraction", false, 0.00f, "minimum allele fraction to be considered in normal, useful for normal sample contaminated with tumor"),
    tumor_f_pretest(easehts::BaseOption::NO_OPTION, "tumor_f_pretest", false, 0.005f, "for computational efficiency, reject sites with allelic fraction below this threshold"),
    min_qscore('q', "min_qscore", false, 5, "threshold for minimum base quality score"),
    gap_events_threshold(easehts::BaseOption::NO_OPTION, "gap_events_threshold", false, 3, "how many gapped events (ins/del) are allowed in proximity to this candidate"),
    heavily_clipped_read_fraction(easehts::BaseOption::NO_OPTION, "heavily_clipped_read_fraction", false, 0.30f, "if this fraction or more of the bases in a read are soft/hard clipped, do not use this read for mutation calling"),
    fraction_mapq0_threshold(easehts::BaseOption::NO_OPTION, "fraction_mapq0_threshold", false, 0.5f, "threshold for determining if there is relatedness between the alt and ref allele read piles"),
    pir_median_threshold(easehts::BaseOption::NO_OPTION, "pir_median_threshold", false, 10, "threshold for clustered read position artifact median"),
    pir_mad_threshold(easehts::BaseOption::NO_OPTION, "pir_mad_threshold", false, 3, "threshold for clustered read position artifact MAD"),
    required_maximum_alt_allele_mapping_quality_score(easehts::BaseOption::NO_OPTION, "required_maximum_alt_allele_mapping_quality_score", false, 20, "required minimum value for tumor alt allele maximum mapping quality score"),
    max_alt_alleles_in_normal_count(easehts::BaseOption::NO_OPTION, "max_alt_alleles_in_normal_count", false, 2, "threshold for maximum alternate allele counts in normal"),
    max_alt_alleles_in_normal_qscore_sum(easehts::BaseOption::NO_OPTION, "max_alt_alleles_in_normal_qscore_sum", false, 20, "threshold for maximum alternate allele quality score sum in normal"),
    max_alt_allele_in_normal_fraction(easehts::BaseOption::NO_OPTION, "max_alt_allele_in_normal_fraction", false, 0.03, "threshold for maximum alternate allele fraction in normal"),
    power_constant_qscore('p', "power_constant_qscore", false, 30, "Phred scale quality score constant to use in power calculations"),
    power_contant_af('P', "power_contant_af", false, 0.3, "Allelic fraction constant to use in power calculations"),
    reference('R', "reference_sequence", true, "Reference sequence file"),
    tumor_files('T', "tumor_file", true, "tumor bam file"),
    normal_files('N', "normal_file", false, "normal bam file"),
    output_file('o', "out", true, "output file"),
    vcf_file('v', "vcf", true, "output vcf file"),
    thread_cnt('t', "nthreads", false, 1, "number of threads"),
    // notice here the default value is false, while in mutect is true
    downsampling(easehts::BaseOption::NO_OPTION, "downsampling", false, "the pileup downsampling")
  {

    easehts::ArgsParser parser;
    parser
      .addOption(enable_qscore_output)
      .addOption(artifact_detection_mode)
      .addOption(force_output)
      .addOption(force_alleles)
      .addOption(only_passing_calls)
      .addOption(initial_tumor_lod_threshold)
      .addOption(tumor_lod_threshold)
      .addOption(fraction_contamination)
      .addOption(minimum_mutation_cell_fraction)
      .addOption(normal_lod_threshold)
      .addOption(normal_artifact_lod_threshold)
      .addOption(strand_artifact_lod_threshold)
      .addOption(strand_artifact_power_threshold)
      .addOption(normal_dbsnp_lod_threshold)
      .addOption(minimum_normal_allele_fraction)
      .addOption(tumor_f_pretest)
      .addOption(min_qscore)
      .addOption(gap_events_threshold)
      .addOption(heavily_clipped_read_fraction)
      .addOption(fraction_mapq0_threshold)
      .addOption(pir_median_threshold)
      .addOption(pir_mad_threshold)
      .addOption(required_maximum_alt_allele_mapping_quality_score)
      .addOption(max_alt_alleles_in_normal_count)
      .addOption(max_alt_alleles_in_normal_qscore_sum)
      .addOption(max_alt_allele_in_normal_fraction)
      .addOption(power_constant_qscore)
      .addOption(power_contant_af)

      .addOption(interval_file)
      .addOption(reference)
      .addOption(tumor_files)
      .addOption(normal_files)
      .addOption(output_file)
      .addOption(vcf_file)
      .addOption(thread_cnt)
      .addOption(downsampling);

    parser.parse(argc, argv);
  }

  easehts::BoolOption enable_qscore_output;
  easehts::BoolOption artifact_detection_mode;
  easehts::BoolOption force_output;
  easehts::BoolOption force_alleles;
  easehts::BoolOption only_passing_calls;
  easehts::FloatOption initial_tumor_lod_threshold;
  easehts::FloatOption tumor_lod_threshold;
  easehts::FloatOption fraction_contamination;
  easehts::FloatOption minimum_mutation_cell_fraction;
  easehts::FloatOption normal_lod_threshold;
  easehts::FloatOption normal_artifact_lod_threshold;
  easehts::FloatOption strand_artifact_lod_threshold;
  easehts::FloatOption strand_artifact_power_threshold;
  easehts::FloatOption normal_dbsnp_lod_threshold;
  easehts::FloatOption minimum_normal_allele_fraction;
  easehts::FloatOption tumor_f_pretest;
  easehts::IntegerOption min_qscore;
  easehts::IntegerOption gap_events_threshold;
  easehts::FloatOption heavily_clipped_read_fraction;
  easehts::FloatOption fraction_mapq0_threshold;
  easehts::DoubleOption pir_median_threshold;
  easehts::DoubleOption pir_mad_threshold;
  easehts::IntegerOption required_maximum_alt_allele_mapping_quality_score;

  /** Parameters for ALT ALLELE IN NORMAL filter **/
  easehts::IntegerOption max_alt_alleles_in_normal_count;
  easehts::IntegerOption max_alt_alleles_in_normal_qscore_sum;
  easehts::DoubleOption max_alt_allele_in_normal_fraction;
  easehts::IntegerOption power_constant_qscore;
  easehts::DoubleOption power_contant_af;

  easehts::StringOption interval_file;
  easehts::StringOption reference;
  easehts::StringListOption tumor_files;
  easehts::StringListOption normal_files;
  easehts::StringOption output_file;
  easehts::StringOption vcf_file;
  easehts::IntegerOption thread_cnt;
  easehts::BoolOption downsampling;
};

} // mutect
} // ncic

#endif
