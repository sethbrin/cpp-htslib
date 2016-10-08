//
// Created by zp on 9/3/16.
//

#ifndef EASEHTSLIB_SAM_BAM_RECORD_H
#define EASEHTSLIB_SAM_BAM_RECORD_H

#include <htslib/sam.h>
#include <string>
#include "utils.h"

namespace ncic {

namespace easehts {

class SAMBAMRecord {
 public:
  /**
   * if is_init is set, init record memory
   * TODO fix design
   */
  SAMBAMRecord(bool is_init=true) {
    if (is_init) {
      raw_record_ = ::bam_init1();
    } else {
      raw_record_ = NULL;
    }
  }

  SAMBAMRecord(SAMBAMRecord&& rhs)
      : raw_record_(rhs.raw_record_) {
    rhs.raw_record_ = NULL;
  }

  bam1_t* GetRawRecord() {
    return raw_record_;
  }

  ~SAMBAMRecord() {
    ::bam_destroy1(raw_record_);
    raw_record_ = NULL;
  }

  SAMBAMRecord(const SAMBAMRecord& rhs) {
    raw_record_ = ::bam_dup1(rhs.raw_record_);
  }

  SAMBAMRecord& operator=(const SAMBAMRecord& rhs) {
    raw_record_ = ::bam_dup1(rhs.raw_record_);
    return *this;
  }

  /*! @function
  @abstract  Get whether the query is on the reverse strand
  @param  b  pointer to an alignment
  @return    boolean true if query is on the reverse strand
  */
  bool IsReverse() {
    return bam_is_rev(raw_record_);
  }

  /*! @function
   @abstract  Get whether the query's mate is on the reverse strand
   @param  b  pointer to an alignment
   @return    boolean true if query's mate on the reverse strand
   */
  bool IsMateReverse() {
    return bam_is_mrev(raw_record_);
  }

  /*! @function
   @abstract  Get the name of the query
   @param  b  pointer to an alignment
   @return    pointer to the name string, null terminated
   */
  const char* GetQueryName() {
    return bam_get_qname(raw_record_);
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

  uint32_t GetSequenceLength() {
    return raw_record_->core.l_qseq;
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

  /*! @function
   @abstract  Get auxiliary data
   @param  b  pointer to an alignment
   @return    pointer to the concatenated auxiliary data
   */
  uint8_t* GetRawAux() {
    return bam_get_aux(raw_record_);
  }

  /*! @function
   @abstract  Get length of auxiliary data
   @param  b  pointer to an alignment
   @return    length of the concatenated auxiliary data
   */
  uint32_t GetRawAuxLength() {
    return bam_get_l_aux(raw_record_);
  }

 private:
  bam1_t* raw_record_;

  // once decode the sequence field, it cached
  std::string cached_sequence_;

};

} // easehts
} // ncic

#endif //EASEHTSLIB_SAM_BAM_RECORD_H
