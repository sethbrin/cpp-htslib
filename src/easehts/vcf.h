//
// Created by zp on 12/15/16.
//

#ifndef EASEHTSLIB_VCF_H_
#define EASEHTSLIB_VCF_H_

#include "easehts/utils.h"

#include <htslib/sam.h>
#include <htslib/vcf.h>

#include <cassert>
#include <memory>

namespace ncic {
namespace easehts {

class VCFHeader {
 public:
  explicit VCFHeader(const char *mode) {
    header_ = ::bcf_hdr_init(mode);
  }

  explicit VCFHeader(bcf_hdr_t* h)
    : header_(h) {}

  ~VCFHeader() {
    ::bcf_hdr_destroy(header_);
  }

  bcf_hdr_t* GetRawHeader() const {
    return header_;
  }

  std::string GetVersion() const {
    return ::bcf_hdr_get_version(header_);
  }

  int GetSamplesCnt() const {
    return bcf_hdr_nsamples(header_);
  }

 private:
  bcf_hdr_t* header_;
};


class VariantContext {
 public:
  VariantContext() {
    record_ = ::bcf_init1();
  }

  ~VariantContext() {
    ::bcf_destroy1(record_);
  }

  bcf1_t* GetRawRecord() const {
    return record_;
  }

 private:
  bcf1_t* record_;
};

class VCFReader {
 public:
  explicit VCFReader(const std::string& filename,
                     const std::string& mode="r") {
    ERROR_COND(!utils::FileExists(filename),
        utils::StringFormatCStr("%s file is not exist!", filename.c_str()));
    filename_ = filename;
    fp_ = ::hts_open(filename.c_str(), mode.c_str());
    header_.reset(new VCFHeader(::bcf_hdr_read(fp_)));
    cur_record_ = nullptr;
  }

  ~VCFReader() {
    int ret;
    if ((ret = ::hts_close(fp_))) {
      fprintf(stderr, "hts_close(%s) non zero status %d\n",
              filename_.c_str(), ret);
      ::exit(ret);
    }
  }

  bool HasNext() {
    VariantContext* record = new VariantContext();

    if (::bcf_read1(fp_, header_->GetRawHeader(),
                    record->GetRawRecord()) >= 0) {
      cur_record_ = record;
      return true;
    }
    delete cur_record_;
    return false;
  }

  VariantContext* Next() {
    return cur_record_;
  }

  const VCFHeader& GetHeader() const {
    return *header_;
  }

 private:
  std::unique_ptr<VCFHeader> header_;
  std::string filename_;
  htsFile* fp_;

  VariantContext* cur_record_;
};

} // easehts
} // ncic

#endif
