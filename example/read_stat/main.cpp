/**
 * example/read_stat/main.cpp
 *
 * stat reads distibution
 */

#include <cstdio>
#include <vector>

#include <easehts/sam_bam_reader.h>
#include <easehts/args_parser.h>

using namespace ncic::easehts;

int main(int argc, char** argv) {
  StringOption input_bam('i', "input", true, "input bam/sam files");
  StringOption output_file('o', "output", false, "output text file");

  ArgsParser parser;
  parser
    .addOption(input_bam)
    .addOption(output_file);
  parser.parse(argc, argv);

  auto work = [&]() {
    FILE* fp = NULL;
    if (output_file.isSet()) {
      fp = fopen(output_file.getValue().c_str(), "w");
    } else {
      fp = stdout;
    }

    SAMBAMTextReader reader(input_bam.getValue());
    SAMBAMRecord record;

    static const int kChunkSize = 100000;
    int start = 0;
    int ref_id = -1;
    int count = 0;
    int pos = 0;
    long long gpos = 0;

    while (reader.HasNext(&record)) {
      int ref = record.GetReferenceId();
      if (ref_id != ref) {
        if (ref_id != -1) {
          fprintf(fp, "%lld,%d\n",
                  gpos, count);
        }
        pos = 0;
        gpos += kChunkSize;
        ref_id = record.GetReferenceId();
        count = 1;
      } else if (record.GetPos() - pos < kChunkSize) {
        count ++;
      } else {
        fprintf(fp, "%lld,%d\n",
                gpos, count);
        count = 1;
        pos += kChunkSize;
        gpos += kChunkSize;
      }
    }
    if (count != 0) {
        fprintf(fp, "%lld,%d\n",
                gpos, count);
    }
    fclose(fp);

  };
  work();
  return 0;
}
