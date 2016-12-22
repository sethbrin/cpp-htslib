#include "easehts/vcf.h"

namespace ncic {
namespace easehts {

void VCFIndexReader::SetRegion(const std::string& region) {
  hts_itr_t *iter = tbx_itr_querys(tbx_idx_, region.c_str());
  if (iter == nullptr) {
    WARN(utils::StringFormatCStr("fail to parse region '%s'", region.c_str()));
    return;
  }
  if (hts_iter_ != nullptr) {
    tbx_itr_destroy(hts_iter_);
  }
  hts_iter_ = iter;
}

void VCFIndexReader::SetRegion(int32_t tid, int32_t begin, int32_t end) {
  hts_itr_t *iter = tbx_itr_queryi(tbx_idx_, tid, begin, end);
  if (iter == nullptr) {
    WARN(utils::StringFormatCStr("fail to parse region '%d:%d-%d'",
                                 tid, begin, end));
    return;
  }
  if (hts_iter_ != nullptr) {
    tbx_itr_destroy(hts_iter_);
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

  kstring_t str = {0,0,0};
  while (tbx_itr_next(fp_, tbx_idx_, hts_iter_, &str) >= 0) {
    if ((vcf_parse(&str, header_->GetRawHeader(),
                   record->GetRawRecord())) >= 0) {
      cur_record_ = record;
      return true;
    } else {
      WARN(utils::StringFormatCStr("Bad format vcf line:%s", str.s));
    }
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

const std::string VCFConstants::kAncestralAlleleKey = "AA";
const std::string VCFConstants::kAlleleCountKey = "AC";
const std::string VCFConstants::kMleAlleleCountKey = "MLEAC";
const std::string VCFConstants::kAlleleFrequenceKey = "AF";
const std::string VCFConstants::kMleAlleleFrequenceKey = "MLEAF";
const std::string VCFConstants::kMlePerSampleAlleleCountKey = "MLPSAC";
const std::string VCFConstants::kMlePerSampleAlleleFractionKey = "MLPAF";
const std::string VCFConstants::kAlleleNumberKey = "AN";
const std::string VCFConstants::kRmsBaseQualityKey = "BQ";
const std::string VCFConstants::kCigarKey = "CIGAR";
const std::string VCFConstants::kDbsnpKey = "DB";
const std::string VCFConstants::kDepthKey = "DP";
const std::string VCFConstants::kDownsampleKey = "DS";
const std::string VCFConstants::kExpectedAlleleCountKey = "EC";
const std::string VCFConstants::kEndKey = "END";
const std::string VCFConstants::kGenotypeFilterKey = "FT";
const std::string VCFConstants::kGenotypeKey = "GT";
const std::string VCFConstants::kGenotypePosteriorsKey = "GP";
const std::string VCFConstants::kGenotypeQualityKey = "GQ";
const std::string VCFConstants::kGenotypeAlleleDepths = "AD";
const std::string VCFConstants::kGenotypePlKey = "PL";


} // easehts
} // ncic
