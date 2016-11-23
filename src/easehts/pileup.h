//
// Created by zp on 11/20/16.
//

#ifndef EASEHTSLIB_PILEUP_H_ 
#define EASEHTSLIB_PILEUP_H_

#include "noncopyable.h"
#include "utils.h"
#include "sam_bam_record.h"

#include <htslib/sam.h>

namespace ncic {
namespace easehts {

/**
 *
 * !@typedef
 * @abstract Structure for one alignment covering the pileup position.
 * @field  b          pointer to the alignment
 * @field  qpos       position of the read base at the pileup site, 0-based
 * @field  indel      indel length; 0 for no indel, positive for ins and negative for del
 * @field  level      the level of the read in the "viewer" mode
 * @field  is_del     1 iff the base on the padded read is a deletion
 * @field  is_head    ???
 * @field  is_tail    ???
 * @field  is_refskip ???
 * @field  aux        ???

 * @discussion See also bam_plbuf_push() and bam_lplbuf_push(). The
 * difference between the two functions is that the former does not
 * set bam_pileup1_t::level, while the later does. Level helps the
 * implementation of alignment viewers, but calculating this has some
 * overhead.
 *
 * typedef struct {
 *     bam1_t *b;
 *     int32_t qpos;
 *     int indel, level;
 *     uint32_t is_del:1, is_head:1, is_tail:1, is_refskip:1, aux:28;
 *     bam_pileup_cd cd; // generic per-struct data, owned by caller.
 * } bam_pileup1_t;
 *
 * XXX PileupElement is just a warper of bam_pileup1_t*,
 * the malloc and free operation is not the business
 */
class PileupElement {
 public:
  explicit PileupElement(const bam_pileup1_t* element)
    : element_(element) {
    qual_ = SAMBAMRecord::GetRawQuality(element_->b)[element_->qpos];
    base_ = utils::SAMUtils::GetSequenceBaseChar(
        SAMBAMRecord::GetRawSequence(element_->b), element_->qpos);
  }

  /*
   * Is this element a deletion w.r.t the reference gnome
   * @return true if this is a deletion, false otherwise*/
  bool IsDeletion() const {
    return element_->is_del;
  }

  uint8_t GetQual() const {
    return qual_;
  }

  bam1_t* GetRead() const {
    return element_->b;
  }

  uint8_t GetBase() const {
    return element_->is_del ? 'D' : base_;
  }

  int GetMappingQuality() const {
    return SAMBAMRecord::GetMapQuality(element_->b);
  }

  int GetOffset() const {
    return element_->qpos;
  }

 public:
  const static char kDeletionBase;
  const static char kDeletionQual;

 private:
  const bam_pileup1_t* element_;

  // TODO cache some field which may compute many times
  uint8_t qual_;
  uint8_t base_;

};

class ReadBackedPileup;
class AbstractReadBackedPileup {
 public:
  AbstractReadBackedPileup(int contig_id, int pos)
    : contig_id_(contig_id),
    pos_(pos) {
    contig_id_pos_ = ((uint64_t)contig_id_ << 32)| pos_;
  }

  AbstractReadBackedPileup()
    : contig_id_(-1),
    pos_(-1),
    contig_id_pos_((uint64_t)-1) {}


  virtual size_t Size() const = 0;
  virtual PileupElement operator[](size_t idx) const = 0;
  virtual PileupElement At(size_t idx) const = 0;

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

  // the pileup site postion
  int contig_id_; // in which contig
  int pos_; // 0-based postion
  uint64_t contig_id_pos_; // combine contig_id and pos

};

/*
 * Encapsulate bam_pileup1_t** structure
 * which contains the pileup reads
 * Each elements is a bam_pileup1_t
 *
 * use the raw pointer to store the pileupelement
 * just for the interface of htslib
 */
class ReadBackedRawPileup : public AbstractReadBackedPileup {
 public:
  ReadBackedRawPileup(const bam_pileup1_t* plp, int size, int contig_id, int pos)
    : AbstractReadBackedPileup(contig_id, pos),
    elements_(plp),
    size_(size) {
  }

  ReadBackedRawPileup()
    : AbstractReadBackedPileup(),
    elements_(nullptr),
    size_(0) {}

  size_t Size() const override {
    return size_;
  }

  PileupElement operator[](size_t idx) const override {
    return PileupElement(elements_ + idx);
  }

  PileupElement At(size_t idx) const override {
    ERROR_COND(idx >= size_,
        utils::StringFormatCStr("out of bound error, the size is %d, and get %d", size_, idx));
    return PileupElement(elements_ + idx);
  }

 private:
  // the PileupElement count
  int size_;

  // the first element position
  const bam_pileup1_t* elements_;
};

/**
 * Use std::vector<PileupElement> to store
 */
class ReadBackedPileup : public AbstractReadBackedPileup {
 public:
  ReadBackedPileup(int contig_id, int pos)
    : AbstractReadBackedPileup(contig_id, pos),
    depth_of_coverage_(kUninitializedCachedIntValue),
    number_of_deletions_(kUninitializedCachedIntValue),
    number_MQ0_reads_(kUninitializedCachedIntValue) {
  }

  ReadBackedPileup()
    : AbstractReadBackedPileup(),
    depth_of_coverage_(kUninitializedCachedIntValue),
    number_of_deletions_(kUninitializedCachedIntValue),
    number_MQ0_reads_(kUninitializedCachedIntValue) {}

  size_t Size() const override {
    return elements_.size();
  }

  PileupElement operator[](size_t idx) const override {
    return elements_[idx];
  }

  PileupElement At(size_t idx) const override {
    return elements_.at(idx);
  }

  void AddElement(PileupElement element) {
    elements_.push_back(element);
  }

  void Clear() {
    elements_.clear();
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
                         std::function<bool (PileupElement element)> pred);

  /**
   * Get the count which satisfy the pred
   */
  int GetPileupByFilterCount(std::function<bool (PileupElement element)> pred) const;

  void GetPositiveStrandPileup(ReadBackedPileup* pPileup);

  void GetNegativeStrandPileup(ReadBackedPileup* pPileup);

  /**
   * Simple useful to count the number of deletion bases in this pileup
   */
  int GetNumberOfDeletions();

  /**
   * Filter the reads with the same read name
   */
  void GetOverlappingFragmentFilteredPileup(ReadBackedPileup* pPileup,
                                            uint8_t ref, bool retain_mismatches);


 private:
  const static int kUninitializedCachedIntValue;

  std::vector<PileupElement> elements_;

  // cache value
  int depth_of_coverage_;
  int number_of_deletions_;
  int number_MQ0_reads_;
};


// TraverseCallback: int(void *data, bam1_t *b)
// status: 0 on success, -1 on end, < -1 on non-recover error
using  TraverseCallback = bam_plp_auto_f;

/*
 * @example
 *
 * while (traverse.HasNext()) {
 *   ReadBackedRawPileup plp = traverse.Next();
 *   printf("contig:%d pos:%d count:%d\n", plp.GetContigId(),
 *          plp.GetPos(), plp.Size());
 * }
 *
 */
class PileupTraverse : public NonCopyable {
 public:
  PileupTraverse(TraverseCallback fun, void* data) {
    iter_ = ::bam_plp_init(fun, data);
  }

  PileupTraverse(PileupTraverse&& rhs)
    : read_backed_pileup_(rhs.read_backed_pileup_),
    iter_(rhs.iter_) {
    rhs.iter_ = nullptr;
  }

  PileupTraverse& operator=(PileupTraverse&& rhs) {
    if (this != &rhs) return *this;

    read_backed_pileup_ = rhs.read_backed_pileup_,
    iter_ = rhs.iter_;
    rhs.iter_ = nullptr;
    return *this;
  }

  ~PileupTraverse() {
    ::bam_plp_destroy(iter_);
  }

  bool HasNext() {
    int contig_id = -1;
    int pos = -1;
    int size = 0;
    const bam_pileup1_t* plp = ::bam_plp_auto(iter_, &contig_id, &pos, &size);

    if (plp != nullptr) {
      read_backed_pileup_ = ReadBackedRawPileup(plp, size, contig_id, pos);
      return true;
    } else {
      read_backed_pileup_ = ReadBackedRawPileup();
      return false;
    }
  }

  const ReadBackedRawPileup& Next() {
    return read_backed_pileup_;
  }

  // the same as Next
  const ReadBackedRawPileup& CurrentPileup() const {
    return read_backed_pileup_;
  }

  void GetNext() {
    HasNext();
  }

 private:
  bam_plp_t iter_; // inner Pileup iterator structure
  ReadBackedRawPileup read_backed_pileup_;
};


} // easehts
} // ncic

#endif
