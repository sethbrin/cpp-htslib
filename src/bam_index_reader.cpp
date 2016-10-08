//
// Created by zp on 9/3/16.
//

#include "bam_index_reader.h"

#include <htslib/sam.h>
#include <htslib/hts.h>

namespace ncic {

namespace easehts {

void BAMIndexReader::AddRegion(const char *region) {
  hts_itr_t *iter = ::sam_itr_querys(index_->GetIndex(),
                                    AbstractReader::header_.GetRawHeader(), region);
  if (iter == NULL) {
    ::fprintf(stderr, "[E::%s] fail to parse region '%s'\n", __func__, region);
    return;
  }

  hts_iters_.push(std::shared_ptr<hts_itr_t>(iter, ::hts_itr_destroy));
}

void BAMIndexReader::AddRegion(int32_t tid, int32_t begin, int32_t end) {
  hts_itr_t *iter = ::sam_itr_queryi(index_->GetIndex(), tid, begin, end);
  if (iter == NULL) {
    ::fprintf(stderr, "[E::%s] fail to load region '%d:%d-%d'\n", __func__, tid, begin, end);
    return;
  }
  hts_iters_.push(std::shared_ptr<hts_itr_t>(iter, ::hts_itr_destroy));
}

bool BAMIndexReader::HasNext(SAMBAMRecord* record) {
  while (!hts_iters_.empty()) {
    std::shared_ptr<hts_itr_t> &front = hts_iters_.front();

    if (sam_itr_next(AbstractReader::fp_, front.get(), record->GetRawRecord()) >= 0) {
      return true;
    }
    // read all the data of front
    hts_iters_.pop();
  }
  return false;
}

} // easehts
} // ncic