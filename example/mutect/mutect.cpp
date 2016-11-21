#include "mutect.h"
#include "locus_read_pile.h"

#include <easehts/genome_loc.h>

#include <atomic>
#include <cmath>
#include <thread>
#include <vector>

namespace ncic {
namespace mutect {

void Worker::Run(const easehts::GenomeLoc& interval) {
  //printf("run: %s\n", interval.ToString().c_str());
  LocusReadPile normal_read_pile(SampleType::NORMAL);
  LocusReadPile tumor_read_pile(SampleType::TUMOR);

  std::vector<easehts::PileupTraverse> normal_traverses;
  std::vector<easehts::PileupTraverse> tumor_traverses;

  int overlap_size = mutect_args_.overlap_size.getValue();
  for (easehts::BAMIndexReader& reader : normal_readers_) {
    InputData input_data;
    reader.SetRegion(interval.GetContigIndex(),
                     interval.GetStart() - overlap_size,
                     interval.GetStop() + overlap_size);
    input_data.reader = &reader;
    normal_traverses.emplace_back(&Worker::Input, &input_data);
  }

  for (easehts::BAMIndexReader& reader : tumor_readers_) {
    InputData input_data;
    reader.SetRegion(interval.GetContigIndex(),
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
    ERROR_COND(cur_contig_id != interval.GetContigIndex(),
               easehts::utils::StringFormatCStr("the contig[%d] is not equal to the interval contig id[%d]",
                                                cur_contig_id, interval.GetContigIndex()));
    if (cur_pos < interval.GetStart() ||
        cur_pos > interval.GetStop()) {
      continue;
    }

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
        tumor_read_pile.Size() == 0) {
      break;
    }

    WARN(easehts::utils::StringFormatCStr("tid:%d pos:%d number of pileup:%d",
                                          cur_contig_id, cur_pos,
                                          normal_read_pile.Size() +
                                          tumor_read_pile.Size()));
    normal_read_pile.Reset();
    tumor_read_pile.Reset();
  }
}

void Worker::PrepareCondidate(int contig_id, int pos,
                              LocusReadPile& normal_read_pile,
                              LocusReadPile& tumor_read_pile) {
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
