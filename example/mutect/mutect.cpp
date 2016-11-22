#include "mutect.h"
#include "locus_read_pile.h"
#include "candidate_mutation.h"

#include <easehts/genome_loc.h>

#include <atomic>
#include <cmath>
#include <thread>
#include <vector>

namespace ncic {
namespace mutect {

const std::string Worker::kValidBases = "ACGT";
const int Worker::kMinQSumQScore = 13;

void Worker::Run(const easehts::GenomeLoc& interval) {
  //printf("run: %s\n", interval.ToString().c_str());

  std::vector<easehts::PileupTraverse> normal_traverses;
  std::vector<easehts::PileupTraverse> tumor_traverses;

  int overlap_size = mutect_args_.overlap_size.getValue();
  for (easehts::BAMIndexReader& reader : normal_readers_) {
    InputData input_data;
    reader.SetRegion(interval.GetContigId(),
                     interval.GetStart() - overlap_size,
                     interval.GetStop() + overlap_size);
    input_data.reader = &reader;
    normal_traverses.emplace_back(&Worker::Input, &input_data);
  }

  for (easehts::BAMIndexReader& reader : tumor_readers_) {
    InputData input_data;
    reader.SetRegion(interval.GetContigId(),
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

    ERROR_COND(cur_contig_id != interval.GetContigId(),
               easehts::utils::StringFormatCStr(
                   "the contig[%d] is not equal to the interval contig id[%d]",
                   cur_contig_id, interval.GetContigId()));

    if (cur_pos < interval.GetStart()) {
      continue;
    }
    if (cur_pos > interval.GetStop()) {
      return;
    }

    easehts::GenomeLoc location(interval.GetContig(), interval.GetContigId(),
                                cur_pos, cur_pos);

    PrepareCondidate(location, min_contig_pos, tumor_traverses, normal_traverses);
  }
}

void Worker::PrepareCondidate(const easehts::GenomeLoc& location,
                              const uint64_t min_contig_pos,
                              const std::vector<easehts::PileupTraverse>& tumor_traverses,
                              const std::vector<easehts::PileupTraverse>& normal_traverses) {
  const char up_ref = std::toupper(reference_.GetSequenceAt(location.GetContig(),
      location.GetStart(), location.GetStop())[0]);
  // only process bases where the reference is [ACGT], because the FASTA
  // for HG18 has N,M and R!
  if (kValidBases.find(up_ref) == std::string::npos) {
    return;
  }
  LocusReadPile normal_read_pile(SampleType::NORMAL, up_ref, mutect_args_.min_qscore.getValue(),
                                 0, true, true, mutect_args_.enable_qscore_output.getValue());
  LocusReadPile tumor_read_pile(SampleType::TUMOR, up_ref, mutect_args_.min_qscore.getValue(),
                                kMinQSumQScore, false, mutect_args_.artifact_detection_mode.getValue(),
                                mutect_args_.enable_qscore_output.getValue());
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
    return;
  }
  WARN(easehts::utils::StringFormatCStr("tid:%d pos:%d number of pileup:%d",
                                        location.GetContigId(), location.GetStart(),
                                        normal_read_pile.Size() +
                                        tumor_read_pile.Size()));

  tumor_read_pile.InitPileups();
  normal_read_pile.InitPileups();
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
