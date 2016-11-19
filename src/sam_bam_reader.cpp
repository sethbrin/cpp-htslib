//
// Created by zp on 11/19/16.
//

#include "sam_bam_reader.h"
#include "utils.h"

#include <string>

namespace ncic {

namespace easehts {

const FileType SRA_TYPE = {"SRA", "sra", ""};
const FileType CRAM_TYPE = {"CRAM", "cram", "crai"};
const FileType BAM_TYPE = {"BAM", "bam", "bai"};
const FileType SAM_TYPE = {"SAM", "sam", ""};


bool SAMBAMTextReader::HasNext(SAMBAMRecord* record) {
  int r = sam_read1(fp_, header_.GetRawHeader(),
                    record->GetRawRecord());
  if (r >= 0) {
    return true;
  } else {
    return false;
  }
}

void BAMIndexReader::SetRegion(const std::string& region) {
  hts_itr_t *iter = ::sam_itr_querys(index_.GetRawIndex(),
                                     header_.GetRawHeader(), region.c_str());
  if (iter == nullptr) {
    WARN(utils::StringFormatCStr("fail to parse region '%s'", region.c_str()));
    return;
  }
  if (hts_iter_ != nullptr) {
    ::hts_itr_destroy(hts_iter_);
  }
  hts_iter_ = iter;
}

void BAMIndexReader::SetRegion(int32_t tid, int32_t begin, int32_t end) {
  hts_itr_t *iter = ::sam_itr_queryi(index_.GetRawIndex(), tid, begin, end);
  if (iter == nullptr) {
    WARN(utils::StringFormatCStr("fail to parse region '%d:%d-%d'", tid, begin, end));
    return;
  }
  if (hts_iter_ != nullptr) {
    ::hts_itr_destroy(hts_iter_);
  }
  hts_iter_ = iter;
}

bool BAMIndexReader::HasNext(SAMBAMRecord* record) {
  if (hts_iter_ == nullptr) {
    // just as SAMBAMTextReader do
    int r = sam_read1(fp_, header_.GetRawHeader(),
                      record->GetRawRecord());
    if (r >= 0) {
      return true;
    } else {
      return false;
    }
  }

  while (::sam_itr_next(fp_, hts_iter_, record->GetRawRecord()) >= 0) {
    return true;
  }
  return false;
}

void BAMIndexBatchReader::AddRegion(const std::string& region) {
  hts_itr_t *iter = ::sam_itr_querys(index_.GetRawIndex(),
                                     header_.GetRawHeader(), region.c_str());
  if (iter == nullptr) {
    WARN(utils::StringFormatCStr("fail to parse region '%s'", region.c_str()));
    return;
  }

  hts_iters_.push(iter);
}

void BAMIndexBatchReader::AddRegion(int32_t tid, int32_t begin, int32_t end) {
  hts_itr_t *iter = ::sam_itr_queryi(index_.GetRawIndex(), tid, begin, end);

  if (iter == nullptr) {
    WARN(utils::StringFormatCStr("fail to parse region '%d:%d-%d'", tid, begin, end));
    return;
  }
  hts_iters_.push(iter);
}

bool BAMIndexBatchReader::HasNext(SAMBAMRecord* record) {
  while (!hts_iters_.empty()) {
    hts_itr_t* front = hts_iters_.front();

    if (sam_itr_next(fp_, front, record->GetRawRecord()) >= 0) {
      return true;
    }
    ::hts_itr_destroy(front);
    // read all the data of front
    hts_iters_.pop();
  }
  return false;
}


} // easehts
} // ncic
