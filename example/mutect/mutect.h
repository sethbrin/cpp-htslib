//
// Created by zp on 11/20/16.
//

#ifndef MUTECT_MUTECT_H_
#define MUTECT_MUTECT_H_

#include "mutect_args.h"
#include "locus_read_pile.h"
#include "power_calculator.h"
#include "candidate_mutation.h"
#include "call_stats_generator.h"

#include <easehts/genome_loc.h>
#include <easehts/noncopyable.h>
#include <easehts/sam_bam_reader.h>
#include <easehts/sam_bam_record.h>
#include <easehts/reference_sequence.h>
#include <easehts/vcf.h>
#include <easehts/gatk/pileup.h>

#include <atomic>
#include <algorithm>
#include <climits>
#include <cfloat>
#include <list>
#include <string>
#include <vector>
#include <memory>

namespace ncic {
namespace mutect {

class Worker : public easehts::NonCopyable {
 public:
  Worker(MutectArgs& mutect_args,
         easehts::IndexedFastaSequenceFile& reference,
         CallStatsGenerator& call_stats_generator)
    : mutect_args_(mutect_args),
    reference_(mutect_args_.reference.getValue()),
    call_stats_generator_(call_stats_generator) {

    // init tumor readers
    std::list<std::string> tumor_values = mutect_args_.tumor_files.getValue();
    for(std::list<std::string>::iterator entry = tumor_values.begin();
         entry != tumor_values.end(); ++entry) {
       tumor_readers_.emplace_back(*entry);
    }

    // init normal readers
    std::list<std::string> normal_values =
      mutect_args_.normal_files.getValue();
    for(std::list<std::string>::iterator entry = normal_values.begin();
        entry != normal_values.end(); ++entry) {
      normal_readers_.emplace_back(*entry);
    }

    // init cosmic traverse
    std::list<std::string> cosmic_values =
      mutect_args_.cosmic_files.getValue();
    for(std::list<std::string>::iterator entry = cosmic_values.begin();
        entry != cosmic_values.end(); ++entry) {
      cosmic_traverses_.emplace_back(
          new easehts::VCFTraverse(new easehts::VCFIndexReader(*entry)));
    }

    // init dbsnp traverse
    std::list<std::string> dbsnp_values =
      mutect_args_.dbsnp_files.getValue();
    for(std::list<std::string>::iterator entry = dbsnp_values.begin();
        entry != dbsnp_values.end(); ++entry) {
      dbsnp_traverses_.emplace_back(
          new easehts::VCFTraverse(new easehts::VCFIndexReader(*entry)));
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

  ~Worker() {
    for (auto& item : cosmic_traverses_) {
      delete item->GetReader();
    }
    for (auto& item : dbsnp_traverses_) {
      delete item->GetReader();
    }
  }

 public:
  const static int kReferenceHalfWindowLength;

 private:
  static int Input(void *data, bam1_t *b);
  // the data
  typedef struct InputData {
    easehts::BAMIndexReader* reader;
  } InputData;

  void PrepareResult(
      const easehts::GenomeLoc& location,
      const uint64_t min_contig_pos,
      const std::vector<easehts::gatk::GATKPileupTraverse>& tumor_traverses,
      const std::vector<easehts::gatk::GATKPileupTraverse>& normal_traverses);

  void PrepareCondidate(
      const char up_ref,
      const easehts::GenomeLoc& location,
      const LocusReadPile& tumor_read_pile,
      const LocusReadPile& normal_read_pile);

  void PrepareVariantContext(const easehts::GenomeLoc& location);
  void PrepareVCFTraverse(const easehts::GenomeLoc& interval);

  void FilterReads(
      const easehts::ReferenceSequence& ref_bases,
      const easehts::GenomeLoc& window,
      const easehts::GenomeLoc& location,
      const easehts::gatk::ReadBackedPileup pileup,
      bool filter_mate_rescue_reads,
      easehts::gatk::ReadBackedPileup* pPileup);

  void PerformRejection(CandidateMutation& candidate);

  const static std::string kValidBases;
  const static int kMinQSumQScore;
  const static int kMaxReadMismatchQualityScoreSum;
  const static char kMappedByMate;

  MutectArgs& mutect_args_;
  easehts::IndexedFastaSequenceFile reference_;
  CallStatsGenerator& call_stats_generator_;

  std::vector<easehts::BAMIndexReader> tumor_readers_;
  std::vector<easehts::BAMIndexReader> normal_readers_;
  bool has_normal_bam_;
  double contaminat_alternate_fraction_;

  std::vector<std::unique_ptr<easehts::VCFTraverse>> cosmic_traverses_;
  std::vector<std::unique_ptr<easehts::VCFTraverse>> dbsnp_traverses_;

  easehts::VariantContext* cosmic_variant_context_;
  easehts::VariantContext* dbsnp_variant_context_;

  TumorPowerCalculator tumor_power_calculator_;
  NormalPowerCalculator normal_novel_site_power_calculator_;
  NormalPowerCalculator normal_db_snp_site_power_calculator_;
  TumorPowerCalculator strand_artifact_power_calculator_;
};

class Mutect : public easehts::NonCopyable {
 public:
  Mutect(int argc, char** argv)
    : mutect_args_(argc, argv),
    reference_(mutect_args_.reference.getValue()),
    call_stats_generator_(mutect_args_.enable_qscore_output.getValue(),
                          mutect_args_.output_file.getValue()) {
  }

  void Run();

 private:
  // XXX the order is very important!!!! you should not change it!!!
  MutectArgs mutect_args_;
  easehts::IndexedFastaSequenceFile reference_;
  CallStatsGenerator call_stats_generator_;
};


/**
 * Collection of Statistical methods and tests used by mutect
 */
class MutectStats : public easehts::NonCopyable {
 public:
  template <typename T>
  static double GetMedian(std::vector<T>& data) {
    std::sort(data.begin(), data.end());
    double result;
    if (data.size() % 2 == 1) {
      // If the number of entries in the lists is not even

      // Get the middle value
      // You must floor the result of the division to drop the remainder
      result = data[data.size() / 2];
    } else {
      result = (data[data.size() / 2] + data[data.size()/2 - 1])/ 2.0;
    }
    return result;
  }

  static double CalculateMAD(const std::vector<int>& data, double median) {
    std::vector<double> dev;
    dev.reserve(data.size());
    for (int d : data) {
      dev.push_back(std::abs(d - median));
    }
    return GetMedian<double>(dev);
  }
};

} // mutect
} // ncic

#endif
