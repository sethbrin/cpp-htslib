#include "parser.h"

#include <easehts/utils.h>

#include <stdio.h>
#include <string>
#include <sys/stat.h>
#include <vector>

int main(int argc, char** argv) {
  StringOption interval('L', "intervals",   true , "interval file");
  StringOption reference('R', "reference_sequence", false, "Reference sequence file");
  StringListOption tumor_files('T', "tumor_file", true, "tumor bam file");
  StringListOption normal_files('N', "normal_file", false, "normal bam file");
  StringOption output_file('o', "out", true, "output file");
  StringOption vcf_file('v', "vcf", true, "output vcf file");

  Parser parser;
  parser.addOption(interval)
    .addOption(reference)
    .addOption(tumor_files)
    .addOption(normal_files)
    .addOption(output_file)
    .addOption(vcf_file);

  std::vector<std::string> otherArguments = parser.parse(argc, argv);


  std::list<std::string> values = tumor_files.getValue();
  for(std::list<std::string>::iterator entry = values.begin();
      entry != values.end();
      ++entry) {
    std::cout << *entry << std::endl;
  }

  return 0;
}
