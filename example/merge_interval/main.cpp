/**
 * example/read_stat/main.cpp
 *
 * stat reads distibution
 */

#include <cstdio>
#include <vector>

#include <easehts/genome_loc.h>
#include <easehts/reference_sequence.h>
#include <easehts/args_parser.h>

using namespace ncic::easehts;

int main(int argc, char** argv) {
  StringOption input_interval('i', "input", true, "input intervals files");
  StringOption output_file('o', "output", false, "output intervals file");
  StringOption reference_file(BaseOption::NO_OPTION, "reference",
                         true, "reference file");

  ArgsParser args_parser;
  args_parser
    .addOption(input_interval)
    .addOption(reference_file)
    .addOption(output_file);
  args_parser.parse(argc, argv);

  IndexedFastaSequenceFile reference(reference_file.getValue());
  auto work = [&]() {
    FILE* fp = NULL;
    if (output_file.isSet()) {
      fp = fopen(output_file.getValue().c_str(), "w");
    } else {
      fp = stdout;
    }

    GenomeLocParser parser(reference.GetSequenceDictionary());
    int padding = 0;
    std::vector<GenomeLoc> intervals =
      IntervalUtils::LoadIntervals(
          parser, input_interval.getValue(), padding);
    for (auto& interval : intervals) {
      fprintf(fp, "%s\t%d\t%d\n", interval.GetContig().c_str(),
              interval.GetStart(), interval.GetStop());
    }
    fclose(fp);
  };

  work();
  return 0;
}
