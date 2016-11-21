#include "mutect.h"

#include <easehts/genome_loc.h>

#include <atomic>
#include <cmath>
#include <thread>
#include <vector>

namespace ncic {
namespace mutect {

void Worker::Run(const easehts::GenomeLoc& interval) {
  printf("run: %s\n", interval.ToString().c_str());
}

void Mutect::Run() {
  easehts::GenomeLocParser parser(reference_.GetSequenceDictionary());
  std::vector<easehts::GenomeLoc> intervals =
    easehts::IntervalUtils::IntervalFileToList(parser,
                                               mutect_args_.interval_file.getValue());

  int thread_cnt = 1;
  if (mutect_args_.thread_cnt.isSet()) {
    thread_cnt = std::max(mutect_args_.thread_cnt.getValue(), thread_cnt);
  }

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
