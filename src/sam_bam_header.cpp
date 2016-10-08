//
// Created by zp on 9/3/16.
//

#include "sam_bam_header.h"

namespace ncic {

namespace easehts {

SAMBAMHeader::SAMBAMHeader(samFile *fp) {
  header_ = ::sam_hdr_read(fp);
}

SAMBAMHeader::~SAMBAMHeader() {
  ::bam_hdr_destroy(header_);
}

} // easehts
} // ncic