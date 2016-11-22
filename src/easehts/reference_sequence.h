//
// Created by zp on 11/18/16.
// Author: Ping Zeng(zengping@ncic.ac.cn)
//

#ifndef EASEHTSLIB_REFERENCE_SEQUENCE_H
#define EASEHTSLIB_REFERENCE_SEQUENCE_H

#include "sam_sequence_dictionary.h"
#include "noncopyable.h"
#include "genome_loc.h"

#include <stdlib.h>
#include <set>
#include <string>
#include <htslib/faidx.h>

namespace ncic {
namespace easehts {

class ReferenceSequence : public NonCopyable {
 public:
  ReferenceSequence()
    : base(NULL),
    len(0) {};

  ReferenceSequence(char* base_, int len_)
    : base(base_),
    len(len_) {}

  ReferenceSequence(ReferenceSequence&& rhs)
    : base(rhs.base),
    len(rhs.len) {
    rhs.base = nullptr;
    len = 0;
  }

  ReferenceSequence& operator=(ReferenceSequence&& rhs) {
    if (this == &rhs) return *this;
    this->base = rhs.base;
    this->len = rhs.len;
    rhs.base = nullptr;
    return *this;
  }

  int Size() const {
    return len;
  }

  char& operator[](size_t idx) {
    return base[idx];
  }

  char& At(size_t idx) {
    ERROR_COND(idx >= len,
        utils::StringFormatCStr("out of bound error, the size is %d, and get %d", len, idx));

    return base[idx];
  }

  ~ReferenceSequence() {
    ::free(base);
  }

  char* base;
  int len;
};

/* Provide core sequence dictionary functionality required by all fasta file readers.
 */
class AbstractFastaSequenceFile : public NonCopyable {
 public:
  AbstractFastaSequenceFile() = default;
  explicit AbstractFastaSequenceFile(std::string file_path);

  AbstractFastaSequenceFile(AbstractFastaSequenceFile&& rhs)
    : file_path_(std::move(rhs.file_path_)),
    ref_header_(std::move(rhs.ref_header_))
  {}

  AbstractFastaSequenceFile& operator=(AbstractFastaSequenceFile&& rhs) {
    file_path_ = std::move(rhs.file_path_);
    ref_header_ = std::move(rhs.ref_header_);
  }

  const SAMSequenceDictionary& GetSequenceDictionary() const {
    return ref_header_;
  }

  std::string GetReferenceFilePath() const {
    return file_path_;
  }

  virtual bool IsIndexed() const {
    return false;
  }

  /**
   * Retrieves the complete sequence described by this contig.
   * @param contig contig whose data should be returned.
   * @return The full sequence associated with this contig.
   */
  virtual ReferenceSequence GetSequence(const std::string& contig) = 0;

  /**
   * Gets the subsequence of the contig in the range [start,stop]
   * @param contig Contig whose subsequence to retrieve.
   * @param start inclusive, 0-based start of region.
   * @param stop inclusive, 0-based stop of region.
   * @return The partial reference sequence associated with this range.
   */
  virtual ReferenceSequence GetSequenceAt(const std::string& contig, int start, int stop) = 0;
  /**
   * Get the dictionary file name of the reference
   * For example: file_path: a.fasta, then the result may be a.dict or
   * a.fasta.dict, if the two files not exist, just return empty string
   * @param reference file name
   * @return dictionary file name
   */
  static std::string FindSequenceDictionaryFileName(const std::string& file_path);
 private:
  // the reference file path
  std::string file_path_;
  // the header of the reference file
  SAMSequenceDictionary ref_header_;

  const static std::set<std::string> kFastaExtensions;
};

/**
 * wrapper of faidx_t struct
 */
class FastaIndex {
 public:
  explicit FastaIndex()
    : fai_(nullptr) {}
  explicit FastaIndex(faidx_t* fai)
    : fai_(fai) {}

  explicit FastaIndex(const char* fn)
    : fai_(LoadFaiIndex(fn)) {}

  ~FastaIndex() {
    DestoryFaiIndex(fai_);
  }

  // You should check if index is null
  bool IsNULL() const {
    return fai_ == nullptr;
  }

  bool HasSequence(const std::string& contig) const {
    CHECK_NOTNULL(fai_);
    return ::faidx_has_seq(fai_, contig.c_str());
  }

  // the count of the contig number
  int Size() const {
    CHECK_NOTNULL(fai_);
    return ::faidx_nseq(fai_);
  }

  // Return name of i-th contig
  const char* GetContigNameAt(size_t idx) const {
    CHECK_NOTNULL(fai_);
    return ::faidx_iseq(fai_, idx);
  }

  int GetContigLength(const std::string& contig) const {
    CHECK_NOTNULL(fai_);
    return ::faidx_seq_len(fai_, contig.c_str());
  }

  /**
   * Fetch the sequence in a region
   * @param contig contig name
   * @param start the start position
   * @param stop the stop position
   * @return ReferenceSequence, ReferenceSequence.len: -2 if seq not present, -1
   * general error, ReferenceSequence.base: Pointer to the sequence; `NULL` on
   * failure
   */
  ReferenceSequence GetSequence(const std::string& contig, int start, int stop) {
    CHECK_NOTNULL(fai_);
    ReferenceSequence ref_seq;
    ref_seq.base = ::faidx_fetch_seq(fai_, contig.c_str(), start, stop, &ref_seq.len);
    return ref_seq;
  }

  /**
   * Fetch the sequence in a region
   * @param region  Region in the format "chr2:20,000-30,000"
   * XXX region is 1-base
   * @return ReferenceSequence, ReferenceSequence.len: -2 if seq not present, -1
   * general error, ReferenceSequence.base: Pointer to the sequence; `NULL` on
   * failure
   *
   */
  ReferenceSequence GetSequence(const std::string& region) {
    CHECK_NOTNULL(fai_);
    ReferenceSequence ref_seq;
    ref_seq.base = ::fai_fetch(fai_, region.c_str(), &ref_seq.len);
    return ref_seq;
  }

  /**
   * Build index for a FASTA or bgzip-compressed FASTA file.
   * @param  fn  FASTA file name
   * @return 0 on success; or -1 on failure
   */
  static int BuildIndex(const char* fn) {
    return ::fai_build(fn);
  }

  /**
   * Load index from "fn.fai".
   *  @param  fn  File name of the FASTA file
   *  @return faidx_t
   */
  static faidx_t* LoadFaiIndex(const char* fn) {
    return ::fai_load(fn);
  }

  /**
   * Destroy a faidx_t struct
   * @param faidx_t struct
   */
  static void DestoryFaiIndex(faidx_t* fai) {
    if (fai == nullptr) return;
    ::fai_destroy(fai);
  }

 private:
  faidx_t* fai_;
};

/**
 * Index FASTA files and extract subsequence.
 *
 * The fai file index columns are:
 *   - chromosome name
 *   - chromosome length: number of bases
 *   - offset: number of bytes to skip to get to the first base
 *     from the beginning of the file, including the length
 *     of the sequence description string (`>chr ..\n`)
 *   - line length: number of bases per line (excluding `\n`)
 *   - binary line length: number of bytes, including `\n`
 */
class IndexedFastaSequenceFile : public AbstractFastaSequenceFile {
 public:
  IndexedFastaSequenceFile() = default;

  explicit IndexedFastaSequenceFile(std::string file_path)
    : AbstractFastaSequenceFile(file_path),
    fasta_index_(file_path.c_str()) {}


  bool IsIndexed() const override {
    return true;
  }

  ReferenceSequence GetSequence(const std::string& contig) override {
    return fasta_index_.GetSequence(contig, 1, fasta_index_.GetContigLength(contig));
  }

  ReferenceSequence GetSequenceAt(const std::string& contig, int start, int stop) override {
    return fasta_index_.GetSequence(contig, start, stop);
  }

 private:
  FastaIndex fasta_index_;

};

} // easehts
} // ncic
#endif
