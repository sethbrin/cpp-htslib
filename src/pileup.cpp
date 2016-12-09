#include "easehts/pileup.h"
#include "easehts/sam_bam_record.h"
#include "easehts/base_utils.h"
#include "easehts/downsampler.h"

#include <assert.h>
#include <string>
#include <unordered_map>
#include <vector>
#include <functional>
#include <memory>

namespace ncic {
namespace easehts {

const char PileupElement::kDeletionBase = 'D';
const char PileupElement::kDeletionQual = (char)16;

const int ReadBackedPileup::kUninitializedCachedIntValue = -1;

void ReadBackedPileup::GetPileupByFilter(ReadBackedPileup* pPileup,
    PileupFilterFun pred) {
  assert(pPileup->Size() == 0);
  int size = elements_.size();
  for (int i = 0; i < size; i++) {
    if (pred(elements_[i])) {
      pPileup->AddElement(elements_[i]);
    }
  }
}

void ReadBackedPileup::GetPileupByAndFilter(ReadBackedPileup* pPileup,
                                            std::vector<PileupFilterFun> filters) {
  assert(pPileup->Size() == 0);
  int size = elements_.size();
  for (int i = 0; i < size; i++) {
    bool flag = true;
    for (const auto& filter : filters) {
      if (!filter(elements_[i])) {
        flag = false;
        break;
      }
    }
    if (flag) {
      pPileup->AddElement(elements_[i]);
    }
  }
}

int ReadBackedPileup::GetPileupByFilterCount(
    PileupFilterFun pred) const {
  int cnt = 0;
  int size = elements_.size();
  for (int i = 0; i < size; i++) {
    if (pred(elements_[i])) {
      cnt ++;
    }
  }
  return cnt;
}

// FIXME here just copy to pPileup and ignore number_of_deletions_
void ReadBackedPileup::GetPileupWithoutDeletions(ReadBackedPileup* pPileup) {

  //auto pred = [](PileupElement element)->bool { return !element.IsDeletion();};
  GetPileupByFilter(pPileup, PileupFilter::IsNotDeletion);
}

/**
 * Get subset of this pileup that contains only bases with
 * quality >= min_base_quality, coming from reads with
 * map quality >= min_map_quality
 * */
void ReadBackedPileup::GetBaseAndMappingFilteredPileup(
    int min_base_quality,
    int min_map_quality,
    ReadBackedPileup* pPileup) {
  //auto pred = [min_base_quality, min_map_quality](PileupElement element)->bool {
    //return SAMBAMRecord::GetMapQuality(element.GetRead()) >= min_map_quality &&
      //(element.IsDeletion() || element.GetQual() >= min_base_quality);
  //};
  GetPileupByFilter(pPileup,
                    std::bind(PileupFilter::IsBaseAndMappingQualityLarge, std::placeholders::_1, min_base_quality, min_map_quality));
}


void ReadBackedPileup::GetBaseFilteredPileup(int min_base_quality,
                                             ReadBackedPileup* pPileup) {
  return GetBaseAndMappingFilteredPileup(min_base_quality, -1, pPileup);
}


int ReadBackedPileup::GetBaseFilteredPileupCount(int min_base_quality) const {
  auto pred = [min_base_quality](PileupElement element)->bool {
    return (element.IsDeletion() || element.GetQual() >= min_base_quality);
  };
  return GetPileupByFilterCount(pred);

}

void ReadBackedPileup::GetMappingFilteredPileup(int min_map_quality,
                                                ReadBackedPileup* pPileup) {
  return GetBaseAndMappingFilteredPileup(-1, min_map_quality, pPileup);
}

void ReadBackedPileup::GetPileupWithoutMappingQualityZeroReads(
    ReadBackedPileup* pPileup) {
  //auto pred = [](PileupElement element)->bool {
    //return SAMBAMRecord::GetMapQuality(element.GetRead()) > 0;
  //};
  GetPileupByFilter(pPileup, PileupFilter::IsMappingQualityLargerThanZero);
}

void ReadBackedPileup::GetPositiveStrandPileup(ReadBackedPileup* pPileup) {
  //auto pred = [](PileupElement element)->bool {
  //  return !SAMBAMRecord::GetReadNegativeStrandFlag(element.GetRead());
  //};
  GetPileupByFilter(pPileup, PileupFilter::IsPositiveStrand);
}

void ReadBackedPileup::GetNegativeStrandPileup(ReadBackedPileup* pPileup) {
  //auto pred = [](PileupElement element)->bool {
  //  return SAMBAMRecord::GetReadNegativeStrandFlag(element.GetRead());
  //};
  GetPileupByFilter(pPileup, PileupFilter::IsNegativeStrand);
}

void ReadBackedPileup::GetOverlappingFragmentFilteredPileup(ReadBackedPileup* pPileup,
                                          uint8_t ref, bool retain_mismatches) {
  int size = elements_.size();
  // store the read_name=>read_idx
  std::unordered_map<std::string, int> filter_map;

  std::vector<bool> elements_to_keep(size, false);

  for (int index=0; index<size; index++) {
    std::string read_name = SAMBAMRecord::GetQueryName(elements_[index].GetRead());

    if (filter_map.find(read_name) == filter_map.end()) {
      filter_map[read_name] = index;
      elements_to_keep[index] = true;
    } else {
      int existing_index = filter_map[read_name];
      const PileupElement& cur_pileup = elements_[index];
      const PileupElement& existing_pileup = elements_[existing_index];

      // if the read disagree at this position
      if (existing_pileup.GetBase() != cur_pileup.GetBase()) {
        if (!retain_mismatches) {
          filter_map.erase(read_name);
          elements_to_keep[existing_index] = false;
        } else {
          // keep the mismatching one
          if (cur_pileup.GetBase() != ref) {
            elements_to_keep[existing_index] = false;
            filter_map[read_name] = index;
            elements_to_keep[index] = true;
          }
        }
      } else {
        // otherwise, keep the element with the higher qualitu score
        if (existing_pileup.GetQual() < cur_pileup.GetQual()) {
          elements_to_keep[existing_index] = false;
          filter_map[read_name] = index;
          elements_to_keep[index] = true;
        }
      }
    }
  }

  // add the element with still keep to pPileup
  for (int index=0; index<size; index++) {
    if (elements_to_keep[index]) {
      pPileup->AddElement(elements_[index]);
    }
  }

}

int ReadBackedPileup::GetNumberofMappingQualityZeroReads() const {
  auto pred = [](PileupElement element)->bool {
    return SAMBAMRecord::GetMapQuality(element.GetRead()) == 0;
  };
  GetPileupByFilterCount(pred);
}

int ReadBackedPileup::GetNumberOfDeletions() {
  if (number_of_deletions_ == kUninitializedCachedIntValue) {
    number_of_deletions_ = 0;
    int size = elements_.size();
    for (int i = 0; i < size; i++) {
      if (elements_[i].IsDeletion()) {
        number_of_deletions_ ++;
      }
    }
  }
  return number_of_deletions_;
}

std::array<int, 4> ReadBackedPileup::GetBaseCounts() const {
  std::array<int, 4> counts{0,0,0,0};

  for (const auto& pile : elements_) {
    if (!pile.IsDeletion()) {
      int index = BaseUtils::SimpleBaseToBaseIndex(pile.GetBase());
      if (index != -1) {
        counts[index]++;
      }
    }
  }
  return counts;
}

int ReadBackedPileup::GetMaxMappingQuals() const {
  int res = 0;
  for (const auto& p : elements_) {
    int qual = SAMBAMRecord::GetMapQuality(p.GetRead());
    if (qual > res) {
      res = qual;
    }
  }
  return res;
}


// GATKPileupTraverse
bool GATKPileupTraverse::HasNext() {
  // 删除buffer_list_中已经超过的出头部分。通过判断尾部是否已经越过
  // 当前coordinate来删除，因为read是不等长的，所以有可能更早的read
  // 没有超出，而晚些的更短read已经超出，
  // TODO: maybe error
  while (!buffer_list_.empty() &&
         (!buffer_list_.front()->IsBeforeEnd(cur_coordianate_))) {
    // free
    //bam_destroy1(buffer_list_.front()->read);
    delete buffer_list_.front();
    buffer_list_.pop_front();
  }

  //if (cur_coordianate_ > 15482579 && cur_coordianate_ < 15482604) {
  //  printf("=====location: %d\n", cur_coordianate_ + 1);
  //  for (auto item : buffer_list_) {
  //    bam1_t* read = item->read;
  //    printf("%s-----%d\n", easehts::SAMBAMRecord::GetQueryName(read),
  //           easehts::SAMBAMRecord::GetSequenceLength(read));
  //  }
  //}

  if (!is_eof_) {
    // 加入新的read
    // 需要判断:
    // 1. read尾在coordinate之前，直接跳过
    // 2. read头在coordinate之前，read尾在coordinate之后，加入buffer_list_
    // 3. read头在coordinate之后， break
    while (true) {
      if (!cur_tracker_->IsBeforeEnd(cur_coordianate_)) {
        if (!GetNextFilteredRead()) {
          break;
        }
        if (cur_tracker_ != nullptr) {
          delete cur_tracker_;
        }
        cur_tracker_ = new PileupTracker(read_);
        continue;
      }
      if (cur_tracker_->IsAfterStart(cur_coordianate_)) {
        // 如果遍历的位点截断了一些read，需要先做一个调整，
        // 让machine移动到当前位点
        int read_start = SAMBAMRecord::GetAlignmentStart(
            cur_tracker_->state_machine.GetRawRead());
        for (int i = 0; i < cur_coordianate_ - read_start; i++) {
          cur_tracker_->StepForwardOnGenome();
        }
        buffer_list_.push_back(cur_tracker_);

        if (!GetNextFilteredRead()) {
          break;
        }

        cur_tracker_ = new PileupTracker(read_);
      } else {
        break;
      }
    }
  }

  // remove the pileup
  for (auto iter = buffer_list_.begin(); iter != buffer_list_.end();) {
    PileupTracker* tracker = *iter;
    // 所以buffer中并不是所有read都可以生成pileup
    // StepForwardOnGenome
    // 走下一格，如果返回false，删除
    if (!tracker->StepForwardOnGenome()) {
      delete tracker;
      iter = buffer_list_.erase(iter);
      continue;
    }
    ++iter;
  }

  if (buffer_list_.size() > downsampler_.GetToCoverage()) {
    downsampler_.DownsampleByAlignmentStart(&buffer_list_);
  }

  // free the pileup element
  FreeReadBackedPileup();
  read_backed_pileup_.reset(new ReadBackedPileup(interval_.GetContigId(),
                                                 cur_coordianate_));
  for (auto iter = buffer_list_.begin(); iter != buffer_list_.end();) {
    PileupTracker* tracker = *iter;
    // LocusIteratorByState.java
    if (SAMBAMRecord::IsBaseInsideAdaptor(
            tracker->read,
            cur_coordianate_)) {
      ++iter;
      continue;
    }


    if (tracker->IsBeforeEnd(cur_coordianate_)) {
      PileupElement element;
      if (tracker->state_machine.MakePileupElement(&element)) {
        read_backed_pileup_->AddElement(element);
      }
    }
    ++iter;
  }

  cur_coordianate_++;
  return true;
}

ReadBackedPileup& GATKPileupTraverse::Next() {
  return *read_backed_pileup_;
}

} // easehts
} // ncic
