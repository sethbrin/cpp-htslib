#include "easehts/pileup.h"
#include "easehts/sam_bam_record.h"
#include "easehts/base_utils.h"

#include <assert.h>
#include <string>
#include <unordered_map>
#include <vector>
#include <functional>

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

} // easehts
} // ncic
