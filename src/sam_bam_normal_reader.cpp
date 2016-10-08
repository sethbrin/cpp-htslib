//
// Created by zp on 8/24/16.
//

#include "sam_bam_normal_reader.h"

#include <htslib/kseq.h>

namespace ncic {

namespace easehts {

void SAMBAMNormalReader::ReadHeader() {
  //SAMTextHeaderCodec codec;
  //file_header_ = codec.Parse(fp_);

}

bool SAMBAMNormalReader::HasNext(SAMBAMRecord* record) {

  int r = sam_read1(AbstractReader::fp_,
                    AbstractReader::header_.GetRawHeader(), record->GetRawRecord());
  if (r >= 0) {
    return true;
  } else {
    return false;
  }

}

} // easehts
} // ncic
