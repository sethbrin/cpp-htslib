//
// Created by zp on 9/3/16.
//

#ifndef EASEHTSLIB_SAM_SEQUENCE_DICTIONARY_H
#define EASEHTSLIB_SAM_SEQUENCE_DICTIONARY_H

#include "utils.h"
#include "noncopyable.h"

#include <htslib/sam.h>
#include <assert.h>
#include <string>
#include <htslib/kseq.h>
#include <htslib/hts.h>
#include <htslib/kstring.h>


namespace ncic {

namespace easehts {

class SAMSequenceDictionary : public NonCopyable {
 public:
  explicit SAMSequenceDictionary() {
    header_ = nullptr;
  }

  explicit SAMSequenceDictionary(samFile *fp) {
    header_ = ::sam_hdr_read(fp);
  }
  explicit SAMSequenceDictionary(bam_hdr_t *header)
    : header_(header) {}

  SAMSequenceDictionary(SAMSequenceDictionary && rhs)
    : header_(rhs.GetRawHeader()) {
    rhs.SetRawHeader(nullptr);
  }

  SAMSequenceDictionary& operator=(SAMSequenceDictionary&& rhs) {
    if (this == &rhs) return *this;
    header_ = rhs.GetRawHeader();
    rhs.SetRawHeader(nullptr);
    return *this;
  }

  ~SAMSequenceDictionary() {
    if (header_ == nullptr) return;
    ::bam_hdr_destroy(header_);
  }

  SAMSequenceDictionary CopyHeader() {
    CHECK_NOTNULL(header_);
    bam_hdr_t* new_header = ::bam_hdr_dup(header_);

    return SAMSequenceDictionary(new_header);
  }

  bam_hdr_t* GetRawHeader() const {
    return header_;
  }

  void SetRawHeader(bam_hdr_t* h) {
    header_ = h;
  }

  int GetSequenceId(const char* ref) const {
    CHECK_NOTNULL(header_);
    return ::bam_name2id(header_, ref);
  }

  int GetSequenceId(const std::string& ref) const {
    CHECK_NOTNULL(header_);
    return ::bam_name2id(header_, ref.c_str());
  }

  bool HasSequence(const char* ref) const {
    CHECK_NOTNULL(header_);
    return ::bam_name2id(header_, ref) != -1;
  }
  bool HasSequence(const std::string& ref) const{
    CHECK_NOTNULL(header_);
    return ::bam_name2id(header_, ref.c_str()) != -1;
  }

  int GetSequenceLen(const char* ref) const {
    CHECK_NOTNULL(header_);
    return header_->target_len[::bam_name2id(header_, ref)];
  }
  int GetSequenceLen(const std::string& ref) const {
    CHECK_NOTNULL(header_);
    return header_->target_len[::bam_name2id(header_, ref.c_str())];
  }

  int WriteToFile(BGZF* fp) {
    CHECK_NOTNULL(header_);
    return ::bam_hdr_write(fp, header_);
  }

  static SAMSequenceDictionary ReadHeader(samFile* fp) {
    kstring_t str;
    bam_hdr_t *h;
    int has_SQ = 0;
    str.l = str.m = 0; str.s = 0;
    while (hts_getline(fp, KS_SEP_LINE, &fp->line) >= 0) {
      if (fp->line.s[0] != '@') break;
      if (fp->line.l > 3 && strncmp(fp->line.s,"@SQ",3) == 0) has_SQ = 1;
      kputsn(fp->line.s, fp->line.l, &str);
      kputc('\n', &str);
    }
    if (str.l == 0) kputsn("", 0, &str);
    h = sam_hdr_parse(str.l, str.s);
    h->l_text = str.l; h->text = str.s;

    return SAMSequenceDictionary(h);
  }

 private:
  bam_hdr_t *header_;
};

} // easehts
} // ncic

#endif //EASEHTSLIB_SAM_SEQUENCE_DICTIONARY_HH
