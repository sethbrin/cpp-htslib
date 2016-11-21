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
    thread_cnt('t', "nthreads", false, "number of threads") {

    easehts::ArgsParser parser;
    parser.addOption(interval_file)
      .addOption(reference)
      .addOption(tumor_files)
      .addOption(normal_files)
      .addOption(output_file)
      .addOption(vcf_file)
      .addOption(thread_cnt);

    parser.parse(argc, argv);
  }

  easehts::StringOption interval_file;
  easehts::StringOption reference;
  easehts::StringListOption tumor_files;
  easehts::StringListOption normal_files;
  easehts::StringOption output_file;
  easehts::StringOption vcf_file;
  easehts::IntegerOption thread_cnt;

};

} // mutect
} // ncic

#endif
