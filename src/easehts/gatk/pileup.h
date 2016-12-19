//
// Created by zp on 12/12/16.
//

#ifndef EASEHTSLIB_PILEUP_GATK_H_
#define EASEHTSLIB_PILEUP_GATK_H_

/**
 * The header file encapsulates the gatk pileup algorithm
 * while pileup.h is a wrapper of samtools pileup
 *
 * The interface of it is the same as pileup.h
 */

#include "easehts/sam_bam_record.h"
#include "easehts/sam_bam_reader.h"
#include "easehts/genome_loc.h"
#include "easehts/utils.h"
#include "easehts/noncopyable.h"

#include "easehts/gatk/pileup_tracker.h"
#include "easehts/gatk/downsampler.h"

#include <memory>
#include <list>

namespace ncic {
namespace easehts {
namespace gatk {

/**
 * GATK pileup element
 */
class PileupElement {
 public:
  PileupElement(SAMBAMRecord* record,
                uint8_t qual, uint8_t base,
                bool is_del, int offset)
    : record_(record),
    qual_(qual),
    base_(base),
    is_del_(is_del),
    offset_(offset) {}

  /*
   * Is this element a deletion w.r.t the reference gnome
   * @return true if this is a deletion, false otherwise*/
  bool IsDeletion() const {
    return is_del_;
  }

  uint8_t GetQual() const {
    return qual_;
  }

  SAMBAMRecord* GetRead() const {
    return record_;
  }

  uint8_t GetBase() const {
    return is_del_ ? 'D' : base_;
  }

  int GetMappingQuality() const {
    return record_->GetMapQuality();
  }

  int GetOffset() const {
    return offset_;
  }

  const static char kDeletionBase;
  const static char kDeletionQual;

 private:
  SAMBAMRecord* record_;
  uint8_t qual_;
  uint8_t base_;
  bool is_del_;
  int offset_;
};

#define GET_PILEUP_BY_FILTER(pileup, pred) \
  assert(pileup->Size() == 0); \
  int size = elements_.size(); \
  for (int i = 0; i < elements_.size(); i++) { \
    if (pred) { \
      pileup->AddElement(elements_[i]); \
    } \
  }

using PileupFilterFun = std::function<bool (PileupElement*)>;
class ReadBackedPileup {
 public:
  ReadBackedPileup(int contig_id, int pos)
    : contig_id_(contig_id),
    pos_(pos),
    contig_id_pos_(((uint64_t)contig_id << 32)| pos),
    depth_of_coverage_(kUninitializedCachedIntValue),
    number_of_deletions_(kUninitializedCachedIntValue),
    number_MQ0_reads_(kUninitializedCachedIntValue) {
    }

  ReadBackedPileup()
    : contig_id_(-1),
    pos_(-1),
    contig_id_pos_(-1),
    depth_of_coverage_(kUninitializedCachedIntValue),
    number_of_deletions_(kUninitializedCachedIntValue),
    number_MQ0_reads_(kUninitializedCachedIntValue) {}

  size_t Size() const {
    return elements_.size();
  }

  PileupElement* operator[](size_t idx) const {
    return elements_[idx];
  }

  PileupElement* At(size_t idx) const {
    return elements_.at(idx);
  }

  void AddElement(PileupElement* element) {
    elements_.push_back(element);
  }

  void Clear() {
    elements_.clear();
  }

  const std::vector<PileupElement*>& GetElements() const {
    return elements_;
  }

  /**
   * Returns a new ReadBackedPileup that is free of deletion spanning reads in
   * this pileup.  Note that this does not copy the data, so both
   * ReadBackedPileups should not be changed.  Doesn't make an unnecessary copy
   * of the pileup (just returns this) if there are no deletions in the pileup.
   */
  void GetPileupWithoutDeletions(ReadBackedPileup* pPileup);

  void GetBaseAndMappingFilteredPileup(int min_base_quality,
                                       int min_map_quality,
                                       ReadBackedPileup* pPileup);


  /**
   * Get subset of this pileup only bases with
   * quality >= min_quality_score
   */
  void GetBaseFilteredPileup(int min_base_quality, ReadBackedPileup* pPileup);

  int GetBaseFilteredPileupCount(int min_base_quality) const;
  /**
   * Get subset of this pileup only bases with
   * mapping quality >= min_quality_score
   */
  void GetMappingFilteredPileup(int min_map_quality, ReadBackedPileup* pPileup);

  /**
   * Get the reads with mapping quality>0
   */
  void GetPileupWithoutMappingQualityZeroReads(ReadBackedPileup* pPileup);

  /**
   * Get the reads with pred true
   */
  void GetPileupByFilter(ReadBackedPileup* pPileup,
                         PileupFilterFun pred);

  /**
   * Get the reads which pass all the filters
   */
  void GetPileupByAndFilter(ReadBackedPileup* pPileup,
                            std::vector<PileupFilterFun> filters);
  /**
   * Get the count which satisfy the pred
   */
  int GetPileupByFilterCount(PileupFilterFun pred) const;

  void GetPositiveStrandPileup(ReadBackedPileup* pPileup);

  void GetNegativeStrandPileup(ReadBackedPileup* pPileup);

  /**
   * Simple useful to count the number of deletion bases in this pileup
   */
  int GetNumberOfDeletions();

  int GetNumberofMappingQualityZeroReads() const;

  /**
   * Filter the reads with the same read name
   */
  void GetOverlappingFragmentFilteredPileup(
      ReadBackedPileup* pPileup,
      uint8_t ref, bool retain_mismatches);

  /**
   * Get the base count
   */
  std::array<int, 4> GetBaseCounts() const ;

  /**
   * Get the max mapping qualities
   */
  int GetMaxMappingQuals() const;

  int GetContigId() const {
    return contig_id_;
  }

  int GetPos() const {
    return pos_;
  }

  // (uint64_t)contig_id_ << 32| pos_;
  uint64_t GetContigPos() const {
    return contig_id_pos_;
  }



 private:
  const static int kUninitializedCachedIntValue;

  std::vector<PileupElement*> elements_;

  // cache value
  int depth_of_coverage_;
  int number_of_deletions_;
  int number_MQ0_reads_;


  // the pileup site postion
  int contig_id_; // in which contig
  int pos_; // 0-based postion
  uint64_t contig_id_pos_; // combine contig_id and pos
};

/**
 * The filters to filter PileupElement
 *
 * when the PileupFilterFun returns true, the element remains,
 * otherwise removed
 */
class PileupFilter : public NonCopyable {
 public:
  PileupFilter() = delete;
  // ignore the element which is deletion
  static bool IsNotDeletion(PileupElement* element) {
    return !element->IsDeletion();
  }

  /**
   * when the quality >= min_base_quality and
   * map quality >= min_map_quality returns true
   */
  static bool IsBaseAndMappingQualityLarge(
      PileupElement* element,
      int min_base_quality, int min_map_quality) {
    bool a =  element->GetRead()->GetMapQuality() >= min_map_quality;
    bool b = (element->IsDeletion() || element->GetQual() >= min_base_quality);
    return element->GetRead()->GetMapQuality() >= min_map_quality &&
      (element->IsDeletion() || element->GetQual() >= min_base_quality);
  }

  static bool IsBaseQualityLarge(
      PileupElement* element, int min_base_quality) {
    return (element->IsDeletion() || element->GetQual() >= min_base_quality);
  }

  static bool IsMappingQualityLarger(
      PileupElement* element,
      int min_map_quality) {
    return element->GetRead()->GetMapQuality() >= min_map_quality;
  }

  static bool IsMappingQualityLargerThanZero(PileupElement* element) {
    return element->GetRead()->GetMapQuality() > 0;
  }


  static bool IsPositiveStrand(PileupElement* element) {
    return !element->GetRead()->GetReadNegativeStrandFlag();
  }

  static bool IsNegativeStrand(PileupElement* element) {
    return element->GetRead()->GetReadNegativeStrandFlag();
  }
};

/**
 * A implementation of GATK Pileup
 */
class GATKPileupTraverse : public NonCopyable {
 public:
  GATKPileupTraverse(BAMIndexReader* reader,
                     const GenomeLoc& traverse_interval,
                     bool is_downsampling=true,
                     int to_coverage = 1000)
  : interval_(traverse_interval),
    reader_(reader) {
    reader_->SetRegion(interval_.GetContigId(), interval_.GetStart(),
                      interval_.GetStop());

    read_backed_pileup_.reset(new ReadBackedPileup());
    cur_coordianate_ = interval_.GetStart();
    if (!GetNextFilteredRead()) {
      is_eof_ = true;
      cur_tracker_ = nullptr;
    } else {
      cur_tracker_ = new PileupTracker(read_);
    }
    if (is_downsampling) {
      downsampler_.reset(new LevelingDownsampler(to_coverage));
    } else {
      downsampler_.reset(nullptr);
    }
  }

  GATKPileupTraverse(GATKPileupTraverse&& rhs) {
    read_backed_pileup_.reset(rhs.read_backed_pileup_.release());
    reader_ = rhs.reader_;
    std::swap(buffer_list_, rhs.buffer_list_);
    std::swap(downsampler_, rhs.downsampler_);
  }

  GATKPileupTraverse& operator=(GATKPileupTraverse&& rhs) {
    if (this != &rhs) return *this;

    read_backed_pileup_.reset(rhs.read_backed_pileup_.release());
    reader_ = rhs.reader_;
    std::swap(buffer_list_, rhs.buffer_list_);
    std::swap(downsampler_, rhs.downsampler_);
    return *this;
  }


  bool HasNext();
  ReadBackedPileup& Next();
  // the same as Next
  const ReadBackedPileup& CurrentPileup() const {
    return *read_backed_pileup_;
  }

  void GetNext() {
    HasNext();
  }

  ~GATKPileupTraverse() {
    int idx = 0;
    while (!buffer_list_.empty()) {
      // free
      delete buffer_list_.front();
      buffer_list_.pop_front();
    }
    FreeReadBackedPileup();
  }



 private:
  void FreeReadBackedPileup() const {
    for (int idx = 0; read_backed_pileup_ &&
         idx < read_backed_pileup_->Size(); idx++) {
      delete (*read_backed_pileup_)[idx];
    }
  }

  bool GetNextFilteredRead() {
    SAMBAMRecord* record = new SAMBAMRecord();
    while (reader_->HasNext(record)) {
      if (record->GetReadUnmappedFlag() ||
          record->GetAlignmentStart() == SAMBAMRecord::NO_ALIGNMENT_START ||
          record->GetNotPrimaryAlignmentFlag() ||
          record->GetDuplicateReadFlag() ||
          record->GetReadFailsVendorQualityCheckFlag()) {
        // read next
      } else {
        // FIXME bug design, should set record to null, otherwise will delete the
        // memory
        read_ = record;
        return true;
      }
    }

    delete record;
    is_eof_ = true;
    return false;
  }


  GenomeLoc interval_;
  BAMIndexReader* reader_;
  SAMBAMRecord* read_;

  // traverse variables
  int cur_coordianate_;
  std::list<PileupTracker*> buffer_list_;
  std::unique_ptr<ReadBackedPileup> read_backed_pileup_;
  PileupTracker* cur_tracker_;
  bool is_eof_ = false;

  std::unique_ptr<LevelingDownsampler> downsampler_;

};




} // gatk
} // easehts
} // ncic

#endif
