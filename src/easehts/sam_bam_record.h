//
// Created by zp on 9/3/16.
//

#ifndef EASEHTSLIB_SAM_BAM_RECORD_H
#define EASEHTSLIB_SAM_BAM_RECORD_H

#include "utils.h"
#include "noncopyable.h"

#include <htslib/sam.h>
#include <climits>
#include <cmath>
#include <functional>
#include <string>

namespace ncic {

namespace easehts {

class CigarElement {
 public:
  CigarElement(char op, int length)
    : operator_(op),
    length_(length) {}

  CigarElement()
  : operator_('D'),
  length_(0) {}

  char GetOperator() const {
    return operator_;
  }

  int GetLength() const {
    return length_;
  }

  std::string ToString() const {
    return std::to_string(length_) + operator_;
  }

  // MIDNSHP=XB
  enum {
    MATCH = 'M',
    INSERTION = 'I',
    DELETION = 'D',
    SKIPPED_REGION = 'N',
    SOFT_CLIP = 'S',
    HARD_CLIP = 'H',
    PADDING = 'P',
    EQUAL = '=',
    MISMATCH = 'X'
  };
 private:
  char operator_;
  int length_;
};

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
    alignment_end_ = kUninitializedCachedIntValue;
    read_name_hashcode_ = 0;
  }

  SAMBAMRecord(SAMBAMRecord&& rhs)
      : raw_record_(rhs.raw_record_) {
    rhs.raw_record_ = nullptr;
    alignment_end_ = kUninitializedCachedIntValue;
    read_name_hashcode_ = 0;
  }

  SAMBAMRecord& operator=(SAMBAMRecord&& rhs) {
    if (this == &rhs) return *this;
    raw_record_ = rhs.raw_record_;
    rhs.raw_record_ = nullptr;
    alignment_end_ = kUninitializedCachedIntValue;
    read_name_hashcode_ = 0;
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
  bool IsReverse() const {
    return bam_is_rev(raw_record_);
  }

  static bool IsReverse(bam1_t* b) {
    return bam_is_rev(b);
  }

  size_t HashCode() const {
    if (read_name_hashcode_ == 0) {
      read_name_hashcode_ = std::hash<std::string>{}(GetQueryName());
    }
    return read_name_hashcode_;
  }

  /*! @function
   @abstract  Get whether the query's mate is on the reverse strand
   @param  b  pointer to an alignment
   @return    boolean true if query's mate on the reverse strand
   */
  bool IsMateReverse() const {
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
  const std::string& GetQueryName() const {
    if (read_name_.empty()) {
      read_name_ = bam_get_qname(raw_record_);
    }
    return read_name_;
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

  int GetInferredInsertSize() const {
    return raw_record_->core.isize;
  }

  static int GetInferredInsertSize(bam1_t* b) {
    return b->core.isize;
  }

  const std::string& GetSequence() {
    if (cached_sequence_.empty()) {
      cached_sequence_.reserve(raw_record_->core.l_qseq);
      const uint8_t* seq = GetRawSequence();
      for (int i = 0; i < raw_record_->core.l_qseq; i++) {
        cached_sequence_.push_back(
            utils::SAMUtils::GetSequenceBaseChar(seq, i));
      }
    }
    return cached_sequence_;
  }

  static std::string GetSequence(bam1_t* b) {
    std::string res;
    res.reserve(b->core.l_qseq);
    const uint8_t* seq = GetRawSequence(b);
    for (int i = 0; i < b->core.l_qseq; i++) {
      res.push_back(utils::SAMUtils::GetSequenceBaseChar(seq, i));
    }
    return res;
  }

  static char GetSequenceAt(bam1_t* b, size_t idx) {
    return utils::SAMUtils::GetSequenceBaseChar(bam_get_seq(b), idx);
  }
  char GetSequenceAt(size_t idx) {
    return utils::SAMUtils::GetSequenceBaseChar(bam_get_seq(raw_record_), idx);
  }

  /*! @function
   @abstract  Get query quality
   @param  b  pointer to an alignment
   @return    pointer to quality string
   */
  uint8_t* GetRawQuality() {
    return bam_get_qual(raw_record_);
  }

  std::string GetQuality() {
    return SAMBAMRecord::GetQuality(raw_record_);
  }

  static uint8_t* GetRawQuality(bam1_t* b) {
    return bam_get_qual(b);
  }

  static std::string GetQuality(bam1_t* b) {
    std::string res;
    uint8_t* qual = GetRawQuality(b);
    for (int i = 0; i < b->core.l_qseq; i++) {
      res.push_back(qual[i]);
    }

    return res;
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
  bool GetReadUnmappedFlag() const {
    return (raw_record_->core.flag & SAMFlag::READ_UNMAPPED) != 0;
  }
  static bool GetReadUnmappedFlag(bam1_t* b) {
    return (b->core.flag & SAMFlag::READ_UNMAPPED) != 0;
  }

  bool GetReadPairedFlag() const {
    return (raw_record_->core.flag & SAMFlag::READ_PAIRED) != 0;
  }
  static bool GetReadPairedFlag(bam1_t* b) {
    return (b->core.flag & SAMFlag::READ_PAIRED) != 0;
  }

  bool GetMateUnmappedFlag() const {
    return (raw_record_->core.flag & SAMFlag::MATE_UNMAPPED) != 0;
  }
  static bool GetMateUnmappedFlag(bam1_t* b) {
    return (b->core.flag & SAMFlag::MATE_UNMAPPED) != 0;
  }

  /**
   * the alignment is not primary (a read having split hits may have multiple
   * primary alignment records).
   */
  bool GetNotPrimaryAlignmentFlag() const {
    return (raw_record_->core.flag & SAMFlag::NOT_PRIMARY_ALIGNMENT) != 0;
  }
  static bool GetNotPrimaryAlignmentFlag(bam1_t* b) {
    return (b->core.flag & SAMFlag::NOT_PRIMARY_ALIGNMENT) != 0;
  }

  /**
   * strand of the query(false for forward; true for reverse strand)
   */
  bool GetReadNegativeStrandFlag() const {
    return (raw_record_->core.flag & SAMFlag::READ_REVERSE_STRAND) != 0;
  }
  static bool GetReadNegativeStrandFlag(bam1_t* b) {
    return (b->core.flag & SAMFlag::READ_REVERSE_STRAND) != 0;
  }

  bool GetMateNegativeStrandFlag() const {
    return (raw_record_->core.flag & SAMFlag::MATE_REVERSE_STRAND) != 0;
  }
  static bool GetMateNegativeStrandFlag(bam1_t* b) {
    return (b->core.flag & SAMFlag::MATE_REVERSE_STRAND) != 0;
  }

  /**
   * the read is either a PCR duplicate or an optical duplicate.
   */
  bool GetDuplicateReadFlag() const {
    return (raw_record_->core.flag & SAMFlag::DUPLICATE_READ) != 0;
  }
  static bool GetDuplicateReadFlag(bam1_t* b) {
    return (b->core.flag & SAMFlag::DUPLICATE_READ) != 0;
  }

  /**
   * the read fails platform/vendor quality checks.
   */
  bool GetReadFailsVendorQualityCheckFlag() const {
    return (raw_record_->core.flag &
            SAMFlag::READ_FAILS_VENDOR_QUALITY_CHECK) != 0;
  }
  static bool GetReadFailsVendorQualityCheckFlag(bam1_t* b) {
    return (b->core.flag & SAMFlag::READ_FAILS_VENDOR_QUALITY_CHECK) != 0;
  }

  /**
   * @return 0-based inclusive leftmost position of the clipped sequence,
   * or 0 if there is no position.
   */
  int GetAlignmentStart() const {
    return raw_record_->core.pos;
  }
  static int GetAlignmentStart(bam1_t* b) {
    return b->core.pos;
  }

  int GetMateAlignmentStart() const {
    return raw_record_->core.mpos;
  }
  static int GetMateAlignmentStart(bam1_t* b) {
    return b->core.mpos;
  }

  int GetAlignmentEnd() const {
    if (alignment_end_ != kUninitializedCachedIntValue) {
      return alignment_end_;
    }
    const std::vector<CigarElement>& cigars = GetCigar();
    alignment_end_ = raw_record_->core.pos - 1;
    for (auto& cigar : cigars) {
      switch(cigar.GetOperator()) {
        case CigarElement::MATCH:
        case CigarElement::DELETION:
        case CigarElement::SKIPPED_REGION:
        case CigarElement::EQUAL:
        case CigarElement::MISMATCH:
          alignment_end_ += cigar.GetLength();
          break;
        default:
          break;
      }
    }
    return alignment_end_;
  }

  static int GetAlignmentEnd(bam1_t* b) {
    std::vector<CigarElement> cigars = ParseRawCigar(b);
    int len = b->core.pos - 1;
    for (auto& cigar : cigars) {
      switch(cigar.GetOperator()) {
        case CigarElement::MATCH:
        case CigarElement::DELETION:
        case CigarElement::SKIPPED_REGION:
        case CigarElement::EQUAL:
        case CigarElement::MISMATCH:
          len += cigar.GetLength();
          break;
        default:
          break;
      }
    }
    return len;
  }


  int GetMapQuality() const {
    return raw_record_->core.qual;
  }
  static int GetMapQuality(bam1_t* b) {
    return b->core.qual;
  }

  uint32_t* GetRawCigar() const {
    return bam_get_cigar(raw_record_);
  }
  static uint32_t* GetRawCigar(bam1_t* b) {
    return bam_get_cigar(b);
  }

  int GetCigarLength() const {
    return raw_record_->core.n_cigar;
  }
  static int GetCigarLength(bam1_t* b) {
    return b->core.n_cigar;
  }

  const std::vector<CigarElement>& GetCigar() const {
    if (cached_cigars_.empty()) {
      cached_cigars_ = ParseRawCigar(raw_record_);
    }
    return cached_cigars_;
  }

  static std::vector<CigarElement> ParseRawCigar(bam1_t* b) {
    std::vector<CigarElement> res;
    uint32_t* cigars = GetRawCigar(b);
    for (int idx=0; idx<b->core.n_cigar; idx++) {
      res.push_back(CigarElement(bam_cigar_opchr(cigars[idx]),
                                 bam_cigar_oplen(cigars[idx])));
    }
    return res;
  }

  std::string GetCigarString() {
    if (cached_cigar_str_.empty()) {
      const std::vector<CigarElement>& elements = GetCigar();
      for (const CigarElement& element : elements) {
        cached_cigar_str_ += element.ToString();
      }
    }
    return cached_cigar_str_;
  }

  static std::string GetCigarString(bam1_t* b) {
    std::vector<CigarElement> elements = ParseRawCigar(b);
    std::string res;
    for (const CigarElement& element : elements) {
      res += element.ToString();
    }
    return res;
  }

  /**
   * Can the adaptor sequence of read be reliably removed from the
   * read base on the alignment of read and its mate?
   *
   * @param read the read to check
   * @return true if it can, false otherwise
   */
  static bool HasWellDefinedFragmentSize(bam1_t* read) {
    if (GetInferredInsertSize(read) == 0) return false;
    if (!GetReadPairedFlag(read)) return false;
    if (GetReadUnmappedFlag(read) || GetMateUnmappedFlag(read)) return false;
    if (GetReadNegativeStrandFlag(read) == GetMateNegativeStrandFlag(read))
      return false;
    if (GetReadNegativeStrandFlag(read)) {
      return GetAlignmentEnd(read) > GetMateAlignmentStart(read);
    } else {
      return GetAlignmentStart(read) <=
        GetMateAlignmentStart(read) + GetInferredInsertSize(read);
    }
  }

  static bool HasWellDefinedFragmentSize(const SAMBAMRecord& read) {
    if (read.GetInferredInsertSize() == 0) return false;
    if (!read.GetReadPairedFlag()) return false;
    if (read.GetReadUnmappedFlag() || read.GetMateUnmappedFlag()) return false;
    if (read.GetReadNegativeStrandFlag() == read.GetMateNegativeStrandFlag())
      return false;
    if (read.GetReadNegativeStrandFlag()) {
      return read.GetAlignmentEnd() > read.GetMateAlignmentStart();
    } else {
      return read.GetAlignmentStart() <=
        read.GetMateAlignmentStart() + read.GetInferredInsertSize();
    }
  }


   /**
    * Finds the adaptor boundary around the read and returns the first
    * base inside the adaptor that is closest to
    * the read boundary. If the read is in the positive strand,
    * this is the first base after the end of the
    * fragment (Picard calls it 'insert'), if the read is in the negative
    * strand, this is the first base before the
    * beginning of the fragment.
    *
    * There are two cases we need to treat here:
    *
    * 1) Our read is in the reverse strand :
    *
    *     <----------------------| *
    *   |--------------------->
    *
    *   in these cases, the adaptor boundary is at the mate start (minus one)
    *
    * 2) Our read is in the forward strand :
    *
    *   |---------------------->   *
    *     <----------------------|
    *
    *   in these cases the adaptor boundary is at the start of the read plus
    *   the inferred insert size (plus one)
    *
    * @param read the read being tested for the adaptor boundary
    * @return the reference coordinate for the adaptor boundary (effectively
    * the first base IN the adaptor, closest to the read.
    * CANNOT_COMPUTE_ADAPTOR_BOUNDARY if the read is unmapped or the
    * MATE_UNMAPPED is mapped to another contig.
    */
  static int GetAdaptorBoundary(bam1_t* read) {
    if (!HasWellDefinedFragmentSize(read)) {
      return CANNOT_COMPUTE_ADAPTOR_BOUNDARY;
    } else if(GetReadNegativeStrandFlag(read)) {
      return GetMateAlignmentStart(read) - 1;
    } else {
      return GetAlignmentStart(read) +
        std::abs(GetInferredInsertSize(read)) + 1;
    }
  }
  static int GetAdaptorBoundary(const SAMBAMRecord& read) {
    if (!HasWellDefinedFragmentSize(read)) {
      return CANNOT_COMPUTE_ADAPTOR_BOUNDARY;
    } else if(read.GetReadNegativeStrandFlag()) {
      return read.GetMateAlignmentStart() - 1;
    } else {
      return read.GetAlignmentStart() +
        std::abs(read.GetInferredInsertSize()) + 1;
    }
  }


  /**
   * is this base inside the adaptor of the read?
   *
   * There are two cases to treat here:
   *
   * 1) Read is in the negative strand => Adaptor boundary is on the left tail
   * 2) Read is in the positive strand => Adaptor boundary is on the right tail
   *
   * Note: We return false to all reads that are UNMAPPED or have an weird
   * big insert size (probably due to mismapping or bigger event)
   *
   * @param read the read to test
   * @param basePos base position in REFERENCE coordinates
   * (not read coordinates)
   * @return whether or not the base is in the adaptor
   */
  static bool IsBaseInsideAdaptor(bam1_t* read, int base_pos) {
    const int adaptor_boundary = GetAdaptorBoundary(read);
    if (adaptor_boundary == CANNOT_COMPUTE_ADAPTOR_BOUNDARY ||
        GetInferredInsertSize(read) > DEFAULT_ADAPTOR_SIZE) {
      return false;
    }
    return GetReadNegativeStrandFlag(read) ?
      base_pos <= adaptor_boundary : base_pos >= adaptor_boundary;
  }
  static bool IsBaseInsideAdaptor(const SAMBAMRecord& read, int base_pos) {
    const int adaptor_boundary = GetAdaptorBoundary(read);
    if (adaptor_boundary == CANNOT_COMPUTE_ADAPTOR_BOUNDARY ||
        read.GetInferredInsertSize() > DEFAULT_ADAPTOR_SIZE) {
      return false;
    }
    return read.GetReadNegativeStrandFlag() ?
      base_pos <= adaptor_boundary : base_pos >= adaptor_boundary;
  }


 public:
  enum {
    CANNOT_COMPUTE_ADAPTOR_BOUNDARY = INT_MIN,
    DEFAULT_ADAPTOR_SIZE = 100,
  };

 private:
  bam1_t* raw_record_;
  static const int kUninitializedCachedIntValue;

  // once decode the sequence field, it cached
  mutable std::string cached_sequence_;
  mutable std::vector<CigarElement> cached_cigars_;
  mutable std::string cached_cigar_str_;
  mutable int alignment_end_;
  mutable std::string read_name_;
  mutable size_t read_name_hashcode_;
};

} // easehts
} // ncic

#endif //EASEHTSLIB_SAM_BAM_RECORD_H
