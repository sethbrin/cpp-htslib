#include <htslib/tbx.h>
#include "easehts/vcf.h"

namespace ncic {
namespace easehts {

void VCFIndexReader::SetRegion(const std::string& region) {
  hts_itr_t *iter = bcf_itr_querys(hts_idx_,
                                   header_->GetRawHeader(), region.c_str());
  if (iter == nullptr) {
    WARN(utils::StringFormatCStr("fail to parse region '%s'", region.c_str()));
    return;
  }
  if (hts_iter_ != nullptr) {
    ::hts_itr_destroy(hts_iter_);
  }
  hts_iter_ = iter;
}

void VCFIndexReader::SetRegion(int32_t tid, int32_t begin, int32_t end) {
  hts_itr_t *iter = bcf_itr_queryi(hts_idx_, tid, begin, end);
  if (iter == nullptr) {
    WARN(utils::StringFormatCStr("fail to parse region '%d:%d-%d'",
                                 tid, begin, end));
    return;
  }
  if (hts_iter_ != nullptr) {
    ::hts_itr_destroy(hts_iter_);
  }
  hts_iter_ = iter;
}

bool VCFIndexReader::HasNext() {
  VariantContext* record = new VariantContext();
  if (hts_iter_ == nullptr) {
    int r = bcf_read1(fp_, header_->GetRawHeader(),
                      record->GetRawRecord());
    if (r >= 0) {
      cur_record_ = record;
      return true;
    }
    delete record;
    return false;
  }

  while (hts_itr_next(hts_get_bgzfp(fp_), hts_iter_, record->GetRawRecord(), 0) >= 0) {
  //while (bcf_itr_next(fp_, hts_iter_, record->GetRawRecord()) >= 0) {
    cur_record_ = record;
    return true;
  }
  delete record;
  return false;

}

void VCFTraverse::SeekFroward(const GenomeLoc& interval) {
  if (interval.GetContigId() == cur_contig_id_) {
    ERROR_COND(interval.GetStart() < cur_pos_,
               utils::StringFormatCStr(
                   "Out of order query: query position %s is located before "
                   "the current position %d:%d",
                   interval.ToString().c_str(),
                   cur_contig_id_, cur_pos_));
    ERROR_COND(interval.GetStop() < cur_query_end_,
               utils::StringFormatCStr(
                   "Unsupported querying sequence: current query interval %s "
                   "ends before the end of previous query interval (%d)",
                   interval.ToString().c_str(), cur_query_end_));
  }
  cur_pos_ = interval.GetStart();
  cur_query_end_ = interval.GetStop();

  if (interval.GetContigId() == cur_contig_id_ &&
      cur_pos_ <= max_pos_) {
    PurgeOutofScopeRecords();
  } else {
    buffer_list_.clear();
    max_pos_ = -1;
    cur_contig_id_ = interval.GetContigId();
  }

  // cur_contig_ and cur_pos_ are set to where we asked to scroll  to
  // if cur_record_ exists, we need not to get next record
  // reader_->Next() equals to cur_record_
  while (cur_record_ || reader_->HasNext()) {
    VariantContext* record = reader_->Next();
    cur_record_ = nullptr;

    GenomeLoc current_contig =
      reader_->GetHeader().CreateOverEntireContig(cur_contig_id_);
    GenomeLoc that_contig = record->GetLocation();

    if (current_contig.IsPast(that_contig)) continue;
    if (current_contig.IsBefore(that_contig)) {
      cur_record_ = record;
      break;
    }

    // we get here if we are on the requested contig
    if (that_contig.GetStop() < cur_pos_) continue;
    if (that_contig.GetStart() > cur_query_end_) {
      cur_record_ = record;
      break;
    }

    if (that_contig.GetStop() > max_pos_) {
      max_pos_ = that_contig.GetStop();
    }

    buffer_list_.push_back(record);
    cur_record_ = nullptr;
  }

}

VariantContext* VCFTraverse::GetFirstRecord() {
  if (!buffer_list_.empty()) {
    return buffer_list_.front();
  }
  return nullptr;
}

void VCFTraverse::GetRecords(std::vector<VariantContext*>* records) {
  records->clear();
  if (!buffer_list_.empty()) {
    for (auto item : buffer_list_) {
      records->push_back(item);
    }
  }
}

void VCFTraverse::PurgeOutofScopeRecords() {
  for (auto iter = buffer_list_.begin(); iter != buffer_list_.end();) {
    if ((*iter)->GetLocation().GetStop() < cur_pos_) {
      delete *iter;
      iter = buffer_list_.erase(iter);
    } else {
      ++iter;
    }
  }
}

} // easehts
} // ncic
