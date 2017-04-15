/**
 * example/read_stat/main.cpp
 *
 * stat reads distibution
 */

#include <cstdio>
#include <vector>
#include <thread>

#include <easehts/sam_bam_reader.h>
#include <easehts/args_parser.h>
#include <easehts/genome_loc.h>
#include <easehts/gatk/pileup.h>
#include <easehts/reference_sequence.h>

using namespace ncic::easehts;

class Worker {
 public:
  Worker(StringOption& reference_file, StringOption& input_bam)
    : input_bam_(input_bam),
    reference_(reference_file.getValue()) {}

  std::pair<uint64_t, uint64_t> Run(const GenomeLoc& interval) {
    BAMIndexReader reader(input_bam_.getValue());
    gatk::GATKPileupTraverse traverse(&reader, interval, false);

    uint64_t count = 0;
    uint64_t total = 0;
    while (traverse.HasNext()) {
      gatk::ReadBackedPileup& pileup = traverse.Next();
      count ++;
      total += pileup.Size();
    }
    return {count, total};
  }
 private:
  StringOption& input_bam_;
  IndexedFastaSequenceFile reference_;
};

int main(int argc, char** argv) {
  StringOption input_bam('i', "input", true, "input bam/sam files");
  StringOption interval_file(BaseOption::NO_OPTION, "intervals",
                             true, "intervals file");
  StringOption reference_file(BaseOption::NO_OPTION, "reference",
                         true, "reference file");

  ArgsParser args_parser;
  args_parser
    .addOption(interval_file)
    .addOption(reference_file)
    .addOption(input_bam);
  args_parser.parse(argc, argv);

  IndexedFastaSequenceFile reference(reference_file.getValue());
  auto work = [&]() {
    GenomeLocParser parser(reference.GetSequenceDictionary());
    std::vector<GenomeLoc> intervals =
      IntervalUtils::LoadIntervals(
          parser, interval_file.getValue(),
          1000);
    int thread_cnt = 48;

    std::atomic<int> interval_index(0);
    std::atomic<int64_t> total(0);
    std::atomic<int64_t> count(0);
    std::vector<std::thread> workers;
    workers.reserve(thread_cnt);
    for (int i = 0; i < thread_cnt; i++) {
      workers.emplace_back([&]() {
         Worker worker(reference_file, input_bam);
         while (interval_index < intervals.size()) {
           size_t index = interval_index++;
           fprintf(stderr, "%d--%d-%d\n", index, intervals[index].GetStart(),
                   intervals[index].GetStop());
           if (index >= intervals.size()) break;

           // As the htslib file reader is [start, end)
           // but the interval is [start, end]
           auto res = worker.Run(GenomeLoc(
                   intervals[index].GetContig(),
                   intervals[index].GetContigId(),
                   intervals[index].GetStart(),
                   intervals[index].GetStop() + 1));
           count += res.first;
           total += res.second;
           fprintf(stderr, "total: %lld count:%lld\n",
                   (int64_t)total, (int64_t)count);
         }
      });
    }

    for (int i = 0; i < thread_cnt; i++) {
      workers[i].join();
    }

    int64_t total_value = total;
    int64_t count_value = count;
    fprintf(stdout, "coverage: %lf\n", (double)total_value/count_value);

  };
  work();
  return 0;
}
