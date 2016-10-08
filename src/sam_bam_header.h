//
// Created by zp on 9/3/16.
//

#ifndef EASEHTSLIB_SAM_BAM_HEADER_H
#define EASEHTSLIB_SAM_BAM_HEADER_H

#include <htslib/sam.h>

namespace ncic {

namespace easehts {

class SAMBAMHeader {
 public:
  explicit SAMBAMHeader(samFile *fp);
  ~SAMBAMHeader();

  bam_hdr_t* GetRawHeader() {
    return header_;
  }

 private:
  bam_hdr_t *header_;
};

} // easehts
} // ncic

#endif //EASEHTSLIB_SAM_BAM_HEADER_H
