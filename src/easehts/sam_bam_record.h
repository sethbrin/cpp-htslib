//
// Created by zp on 9/3/16.
//

#ifndef EASEHTSLIB_SAM_BAM_RECORD_H
#define EASEHTSLIB_SAM_BAM_RECORD_H

#include "utils.h"
#include "noncopyable.h"

#include <htslib/sam.h>
#include <string>

namespace ncic {

namespace easehts {

enum SAMFlag {
  // "Template having multiple segments in sequencing"
  READ_PAIRED = 0x1,
  // Each segment properly aligned according to the aligner
  PROPER_PAIR = 0x2,
  // Segment unmapped
  READ_UNMAPPED = 0x4,
  // Next segment in the template unmapped
  MATE_UNMAPPED = 0x8,
  // SEQ being reverse complemented
  READ_REVERSE_STRAND = 0x10,
  // SEQ of the next segment in the template being reverse complemented
  MATE_REVERSE_STRAND = 0x20,
  // The first segment in the template
  FIRST_OF_PAIR = 0x40,
  // The last segment in the template
  SECOND_OF_PAIR = 0x80,
  // Secondary alignment
  NOT_PRIMARY_ALIGNMENT = 0x100,
  // Not passing quality controls
  READ_FAILS_VENDOR_QUALITY_CHECK = 0x200,
  // PCR or optical duplicate
  DUPLICATE_READ = 0x400,
  // Supplementary alignment
  SUPPLEMENTARY_ALIGNMENT = 0x800
};

class SAMBAMRecord : public NonCopyable {
 public:
  enum {
    NO_ALIGNMENT_START = -1
  };
  /**
   * if is_init is set, init record memory
   * TODO fix design
   */
  SAMBAMRecord(bool is_init=true) {
    if (is_init) {
      raw_record_ = ::bam_init1();
    } else {
      raw_record_ = nullptr;
    }
  }

  SAMBAMRecord(SAMBAMRecord&& rhs)
      : raw_record_(rhs.raw_record_) {
    rhs.raw_record_ = nullptr;
  }

  SAMBAMRecord& operator=(SAMBAMRecord&& rhs) {
    if (this == &rhs) return *this;
    raw_record_ = rhs.raw_record_;
    rhs.raw_record_ = nullptr;
    return *this;
  }

  bam1_t* GetRawRecord() {
    return raw_record_;
  }

  void SetRawRecord(bam1_t* b) {
    raw_record_ = b;
  }

  ~SAMBAMRecord() {
    ::bam_destroy1(raw_record_);
    raw_record_ = nullptr;
  }

  SAMBAMRecord Copy() {
    SAMBAMRecord record;
    record.raw_record_ = ::bam_dup1(raw_record_);
    return record;
  }

  /*! @function
  @abstract  Get whether the query is on the reverse strand
  @param  b  pointer to an alignment
  @return    boolean true if query is on the reverse strand
  */
  bool IsReverse() {
    return bam_is_rev(raw_record_);
  }

  static bool IsReverse(bam1_t* b) {
    return bam_is_rev(b);
  }

  /*! @function
   @abstract  Get whether the query's mate is on the reverse strand
   @param  b  pointer to an alignment
   @return    boolean true if query's mate on the reverse strand
   */
  bool IsMateReverse() {
    return bam_is_mrev(raw_record_);
  }

  static bool IsMateReverse(bam1_t* b) {
    return bam_is_mrev(b);
  }

  /*! @function
   @abstract  Get the name of the query
   @param  b  pointer to an alignment
   @return    pointer to the name string, null terminated
   */
  const char* GetQueryName() {
    return bam_get_qname(raw_record_);
  }

  static const char* GetQueryName(bam1_t* b) {
    return bam_get_qname(b);
  }

  /*! @function
   @abstract  Get query sequence
   @param  b  pointer to an alignment
   @return    pointer to sequence

   @discussion Each base is encoded in 4 bits: 1 for A, 2 for C, 4 for G,
   8 for T and 15 for N. Two bases are packed in one byte with the base
   at the higher 4 bits having smaller coordinate on the read. It is
   recommended to use bam_seqi() macro to get the base.
   */
  uint8_t* GetRawSequence() {
    return bam_get_seq(raw_record_);
  }

  static uint8_t* GetRawSequence(bam1_t* b) {
    return bam_get_seq(b);
  }

  uint32_t GetSequenceLength() {
    return raw_record_->core.l_qseq;
  }

  static uint32_t GetSequenceLength(bam1_t* b) {
    return b->core.l_qseq;
  }

  std::string GetSequence() {
    if (cached_sequence_.empty()) {
      cached_sequence_.reserve(raw_record_->core.l_qseq);
      const uint8_t* seq = GetRawSequence();
      for (int i = 0; i < raw_record_->core.l_qseq; i++) {
        cached_sequence_.push_back(utils::SAMUtils::GetSequenceBaseChar(seq, i));
      }
    }
    return cached_sequence_;
  }

  /*! @function
   @abstract  Get query quality
   @param  b  pointer to an alignment
   @return    pointer to quality string
   */
  uint8_t* GetRawQuality() {
    return bam_get_qual(raw_record_);
  }

  static uint8_t* GetRawQuality(bam1_t* b) {
    return bam_get_qual(b);
  }

  /*! @function
   @abstract  Get auxiliary data
   @param  b  pointer to an alignment
   @return    pointer to the concatenated auxiliary data
   */
  uint8_t* GetRawAux() {
    return bam_get_aux(raw_record_);
  }

  static uint8_t* GetRawAux(bam1_t* b) {
    return bam_get_aux(b);
  }

  /*! @function
   @abstract  Get length of auxiliary data
   @param  b  pointer to an alignment
   @return    length of the concatenated auxiliary data
   */
  uint32_t GetRawAuxLength() {
    return bam_get_l_aux(raw_record_);
  }

  static uint32_t GetRawAuxLength(bam1_t* b) {
    return bam_get_l_aux(b);
  }

  /**
   * the query sequence itself is unmapped.
   */
  bool GetReadUnmappedFlag() {
    return (raw_record_->core.flag & SAMFlag::READ_UNMAPPED) != 0;
  }
  static bool GetReadUnmappedFlag(bam1_t* b) {
    return (b->core.flag & SAMFlag::READ_UNMAPPED) != 0;
  }

  /**
   * the alignment is not primary (a read having split hits may have multiple primary alignment records).
   */
  bool GetNotPrimaryAlignmentFlag() {
    return (raw_record_->core.flag & SAMFlag::NOT_PRIMARY_ALIGNMENT) != 0;
  }
  static bool GetNotPrimaryAlignmentFlag(bam1_t* b) {
    return (b->core.flag & SAMFlag::NOT_PRIMARY_ALIGNMENT) != 0;
  }

  /**
   * strand of the query(false for forward; true for reverse strand)
   */
  bool GetReadNegativeStrandFlag() {
    return (raw_record_->core.flag & SAMFlag::READ_REVERSE_STRAND) != 0;
  }
  static bool GetReadNegativeStrandFlag(bam1_t* b) {
    return (b->core.flag & SAMFlag::READ_REVERSE_STRAND) != 0;
  }

  /**
   * the read is either a PCR duplicate or an optical duplicate.
   */
  bool GetDuplicateReadFlag() {
    return (raw_record_->core.flag & SAMFlag::DUPLICATE_READ) != 0;
  }
  static bool GetDuplicateReadFlag(bam1_t* b) {
    return (b->core.flag & SAMFlag::DUPLICATE_READ) != 0;
  }

  /**
   * the read fails platform/vendor quality checks.
   */
  bool GetReadFailsVendorQualityCheckFlag() {
    return (raw_record_->core.flag & SAMFlag::READ_FAILS_VENDOR_QUALITY_CHECK) != 0;
  }
  static bool GetReadFailsVendorQualityCheckFlag(bam1_t* b) {
    return (b->core.flag & SAMFlag::READ_FAILS_VENDOR_QUALITY_CHECK) != 0;
  }

  /**
   * @return 0-based inclusive leftmost position of the clipped sequence, or 0 if there is no position.
   */
  int GetAlignmentStart() {
    return raw_record_->core.pos;
  }
  static int GetAlignmentStart(bam1_t* b) {
    return b->core.pos;
  }

  int GetMapQuality() {
    return raw_record_->core.qual;
  }
  static int GetMapQuality(bam1_t* b) {
    return b->core.qual;
  }
 private:
  bam1_t* raw_record_;

  // once decode the sequence field, it cached
  std::string cached_sequence_;

};

} // easehts
} // ncic

#endif //EASEHTSLIB_SAM_BAM_RECORD_H
