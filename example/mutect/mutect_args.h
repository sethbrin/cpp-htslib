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
    reference('R', "reference_sequence", true, "Reference sequence file"),
    tumor_files('T', "tumor_file", true, "tumor bam file"),
    normal_files('N', "normal_file", false, "normal bam file"),
    output_file('o', "out", true, "output file"),
    vcf_file('v', "vcf", true, "output vcf file"),
    thread_cnt('t', "nthreads", false, 1, "number of threads"),
    overlap_size('O', "overlap_size", false, 3000, "the size of overlap_size, default is 3000"),
    min_qscore('q', "min_qscore", false, 5, "threshold for minimum base quality score"),
    artifact_detection_mode('a', "artifact_detection_mode", false, "used when running the caller on a normal (as if it were a tumor) to detect artifacts"),
    enable_qscore_output('e', "enable_qscore_output", false, "output quality scores of all bases used for ref and alt")
  {

    easehts::ArgsParser parser;
    parser.addOption(interval_file)
      .addOption(reference)
      .addOption(tumor_files)
      .addOption(normal_files)
      .addOption(output_file)
      .addOption(vcf_file)
      .addOption(thread_cnt)
      .addOption(overlap_size)
      .addOption(min_qscore)
      .addOption(artifact_detection_mode)
      .addOption(enable_qscore_output);

    parser.parse(argc, argv);
  }

  easehts::StringOption interval_file;
  easehts::StringOption reference;
  easehts::StringListOption tumor_files;
  easehts::StringListOption normal_files;
  easehts::StringOption output_file;
  easehts::StringOption vcf_file;
  easehts::IntegerOption thread_cnt;
  easehts::IntegerOption overlap_size;
  easehts::IntegerOption min_qscore;
  easehts::BoolOption artifact_detection_mode;
  easehts::BoolOption enable_qscore_output;
};

} // mutect
} // ncic

#endif
