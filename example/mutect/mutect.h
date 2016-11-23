//
// Created by zp on 11/20/16.
//

#ifndef MUTECT_MUTECT_H_
#define MUTECT_MUTECT_H_

#include "mutect_args.h"
#include "locus_read_pile.h"

#include <easehts/genome_loc.h>
#include <easehts/noncopyable.h>
#include <easehts/sam_bam_reader.h>
#include <easehts/sam_bam_record.h>
#include <easehts/reference_sequence.h>

#include <list>
#include <string>
#include <vector>
#include <atomic>

namespace ncic {
namespace mutect {

class Worker : public easehts::NonCopyable {
 public:
  Worker(MutectArgs& mutect_args,
         easehts::IndexedFastaSequenceFile& reference)
    : mutect_args_(mutect_args),
    reference_(reference) {

    // init tumor readers
    std::list<std::string> tumor_values = mutect_args_.tumor_files.getValue();
    for(std::list<std::string>::iterator entry = tumor_values.begin();
         entry != tumor_values.end(); ++entry) {
       tumor_readers_.emplace_back(*entry);
    }

    // init normal readers
    std::list<std::string> normal_values = mutect_args_.normal_files.getValue();
    for(std::list<std::string>::iterator entry = normal_values.begin();
        entry != normal_values.end(); ++entry) {
      normal_readers_.emplace_back(*entry);
    }

    has_normal_bam_ = normal_readers_.size() != 0;
  }

  void Run(const easehts::GenomeLoc& interval);

 private:
  static int Input(void *data, bam1_t *b);
  // the data
  typedef struct InputData {
    easehts::BAMIndexReader* reader;
  } InputData;

  void PrepareResult(const easehts::GenomeLoc& location,
                     const uint64_t min_contig_pos,
                     const std::vector<easehts::PileupTraverse>& tumor_traverses,
                     const std::vector<easehts::PileupTraverse>& normal_traverses);

  void PrepareCondidate(const LocusReadPile& tumor_read_pile,
                        const LocusReadPile& normal_read_pile);


  const static std::string kValidBases;
  const static int kMinQSumQScore;

  MutectArgs& mutect_args_;
  easehts::IndexedFastaSequenceFile& reference_;

  std::vector<easehts::BAMIndexReader> tumor_readers_;
  std::vector<easehts::BAMIndexReader> normal_readers_;
  bool has_normal_bam_;
};

class Mutect : public easehts::NonCopyable {
 public:
  Mutect(int argc, char** argv)
    : mutect_args_(argc, argv),
    reference_(mutect_args_.reference.getValue()) {
  }

  void Run();

 private:
  // XXX the order is very important!!!! you should not change it!!!
  MutectArgs mutect_args_;
  easehts::IndexedFastaSequenceFile reference_;
};

} // mutect
} // ncic

#endif
