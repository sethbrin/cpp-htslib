//
// Created by zp on 11/20/16.
//

#include "noncopyable.h"
#include "utils.h"

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
    : element_(element) {}

 private:
  const bam_pileup1_t* element_;
};

/*
 * Encapsulate bam_pileup1_t** structure
 * which contains the pileup reads
 * Each elements is a bam_pileup1_t
 */
class ReadBackedPileup {
 public:
  ReadBackedPileup(const bam_pileup1_t* plp, int size, int contig_id, int pos)
    : elements_(plp),
    size_(size),
    contig_id_(contig_id),
    pos_(pos) {
    contig_id_pos_ = ((uint64_t)contig_id_ << 32)| pos_;
  }

  ReadBackedPileup()
    : elements_(nullptr),
    size_(0),
    contig_id_(-1),
    pos_(-1),
    contig_id_pos_((uint64_t)-1) {}

  size_t Size() const {
    return size_;
  }

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

  PileupElement operator[](size_t idx) const {
    return PileupElement(elements_ + idx);
  }

  PileupElement At(size_t idx) const {
    ERROR_COND(idx >= size_,
        utils::StringFormatCStr("out of bound error, the size is %d, and get %d", size_, idx));
    return PileupElement(elements_ + idx);
  }

 private:
  // the first element position
  const bam_pileup1_t* elements_;
  // the PileupElement count
  int size_;

  // the pileup site postion
  int contig_id_; // in which contig
  int pos_; // 0-based postion
  uint64_t contig_id_pos_; // combine contig_id and pos
};

// TraverseCallback: int(void *data, bam1_t *b)
// status: 0 on success, -1 on end, < -1 on non-recover error
using  TraverseCallback = bam_plp_auto_f;

/*
 * @example
 *
 * while (traverse.HasNext()) {
 *   ReadBackedPileup plp = traverse.Next();
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
      read_backed_pileup_ = ReadBackedPileup(plp, size, contig_id, pos);
      return true;
    } else {
      read_backed_pileup_ = ReadBackedPileup();
      return false;
    }
  }

  const ReadBackedPileup& Next() {
    return read_backed_pileup_;
  }

  // the same as Next
  const ReadBackedPileup& CurrentPileup() {
    return read_backed_pileup_;
  }

  void GetNext() {
    HasNext();
  }

 private:
  bam_plp_t iter_; // inner Pileup iterator structure
  ReadBackedPileup read_backed_pileup_;
};


} // easehts
} // ncic
