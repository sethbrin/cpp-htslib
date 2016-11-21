//
// Created by zp on 11/20/16.
//

#ifndef MUTECT_MUTECT_H_
#define MUTECT_MUTECT_H_

#include "mutect_args.h"

#include <easehts/genome_loc.h>
#include <easehts/noncopyable.h>
#include <easehts/reference_sequence.h>
#include <easehts/sam_bam_reader.h>

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
  }

  void Run(const easehts::GenomeLoc& interval);

 private:
  MutectArgs& mutect_args_;
  easehts::IndexedFastaSequenceFile& reference_;

  std::vector<easehts::BAMIndexReader> tumor_readers_;
  std::vector<easehts::BAMIndexReader> normal_readers_;
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
