//
// Created by zp on 11/20/16.
//

#ifndef MUTECT_MUTECT_H_
#define MUTECT_MUTECT_H_

#include "mutect_args.h"
#include "locus_read_pile.h"
#include "power_calculator.h"

#include <easehts/genome_loc.h>
#include <easehts/noncopyable.h>
#include <easehts/sam_bam_reader.h>
#include <easehts/sam_bam_record.h>
#include <easehts/reference_sequence.h>

#include <climits>
#include <cfloat>
#include <atomic>
#include <list>
#include <string>
#include <vector>

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
    if (!has_normal_bam_) {
      mutect_args_.normal_lod_threshold.setValue(-1 * FLT_MAX);
      mutect_args_.normal_dbsnp_lod_threshold.setValue(-1 * FLT_MAX);
      mutect_args_.normal_artifact_lod_threshold.setValue(FLT_MAX);
    }

    contaminat_alternate_fraction_ = std::max(
        mutect_args_.minimum_mutation_cell_fraction.getValue(),
        mutect_args_.fraction_contamination.getValue());

    // coverage related initialization
    double power_constant_eps = std::pow(10,
        -1 * (mutect_args_.power_constant_qscore.getValue()/10));

    tumor_power_calculator_ = TumorPowerCalculator(
        power_constant_eps,
        mutect_args_.tumor_lod_threshold.getValue(),
        contaminat_alternate_fraction_);
    normal_novel_site_power_calculator_ = NormalPowerCalculator(
        power_constant_eps,
        mutect_args_.normal_lod_threshold.getValue());
    normal_db_snp_site_power_calculator_ = NormalPowerCalculator(
        power_constant_eps,
        mutect_args_.normal_dbsnp_lod_threshold.getValue());
    strand_artifact_power_calculator_ = TumorPowerCalculator(
        power_constant_eps,
        mutect_args_.strand_artifact_lod_threshold.getValue(), 0.0f);

    // to force output, all we have to do is lower the initial tumor lod threshold
    if (mutect_args_.force_output.getValue()) {
      mutect_args_.initial_tumor_lod_threshold.setValue(-FLT_MAX);
    }

    // initialize the vcf output

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

  void PrepareCondidate(
      const char up_ref,
      const easehts::GenomeLoc& location,
      const LocusReadPile& tumor_read_pile,
      const LocusReadPile& normal_read_pile);


  const static std::string kValidBases;
  const static int kMinQSumQScore;

  MutectArgs& mutect_args_;
  easehts::IndexedFastaSequenceFile& reference_;

  std::vector<easehts::BAMIndexReader> tumor_readers_;
  std::vector<easehts::BAMIndexReader> normal_readers_;
  bool has_normal_bam_;
  double contaminat_alternate_fraction_;

  TumorPowerCalculator tumor_power_calculator_;
  NormalPowerCalculator normal_novel_site_power_calculator_;
  NormalPowerCalculator normal_db_snp_site_power_calculator_;
  TumorPowerCalculator strand_artifact_power_calculator_;
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
