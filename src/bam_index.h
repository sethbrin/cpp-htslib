//
// Created by zp on 9/3/16.
//

#ifndef EASEHTSLIB_BAM_INDEX_H
#define EASEHTSLIB_BAM_INDEX_H

#include "noncopyable.h"

#include <string>
#include <htslib/hts.h>

namespace ncic {

namespace easehts {

class BAMIndex : public NonCopyable {
 public:
  BAMIndex(const std::string& filename)
      : filename_(filename) {
    hts_idx_ = ::hts_idx_load(filename_.c_str(), HTS_FMT_BAI);
  }

  BAMIndex(const std::string& filename, const std::string& ext)
      : filename_(filename),
        filename_ext_(ext) {
    hts_idx_ = ::hts_idx_load2(filename_.c_str(), filename_ext_.c_str());
  }

  std::string GetFileName() {
    return filename_;
  }

  std::string GetFileExtName() {
    return filename_ext_;
  }

  hts_idx_t* GetIndex() {
    return hts_idx_;
  }

  ~BAMIndex() {
    ::hts_idx_destroy(hts_idx_);
  }

 private:
  std::string filename_; // the data filename
  std::string filename_ext_; // the default is .bai/.csi
  hts_idx_t* hts_idx_;
};

} // easehts
} // ncic

#endif //EASEHTSLIB_BAM_INDEX_H
