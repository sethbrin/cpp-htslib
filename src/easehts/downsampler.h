//
// Created by zp on 12/8/16.
//

#ifndef EASEHTSLIB_DOWNSAMPLER_H_
#define EASEHTSLIB_DOWNSAMPLER_H_

#include "utils.h"
#include "pileup_tracker.h"

#include <list>
#include <vector>

namespace ncic {
namespace easehts {

class LevelingDownsampler {
 public:
  /**
   * Construct a LevelingDownsampler
   *
   * @param target_size_ the sum of the sizes of all individual Lists this downsampler is fed may not exceed
   * this value -- if it does, items are removed from Lists evenly until the total size is <= this value
   *
   * @param min_element_per_stack_ no stack will be reduced below this size during downsampling.  That is,
   * if a stack has only 3 elements and minElementsPerStack is 3, no matter what
   * we'll not reduce this stack below 3.
   */
  LevelingDownsampler(int target_size, int min_element_per_stack) {
    ERROR_COND(target_size < 0,
               utils::StringFormatCStr(
                   "target_size must be >=0 but got %d", target_size));
    ERROR_COND(min_element_per_stack < 0,
               utils::StringFormatCStr(
                   "target_size must be >=0 but got %d", min_element_per_stack));
    target_size_ = target_size;
    min_element_per_stack_ = min_element_per_stack;
  }

  explicit LevelingDownsampler(int target_size)
    : LevelingDownsampler(target_size, 1) {}

  LevelingDownsampler()
    : LevelingDownsampler(0, 1) {}

  LevelingDownsampler(LevelingDownsampler&& rhs)
    : target_size_(rhs.target_size_),
    min_element_per_stack_(rhs.min_element_per_stack_) {}

  LevelingDownsampler& operator=(LevelingDownsampler&& rhs) {
    if (this == &rhs) return *this;

    target_size_ = rhs.target_size_,
    min_element_per_stack_ = rhs.min_element_per_stack_;
    return *this;
  }

  int DownsampleByAlignmentStart(std::list<PileupTracker*>* tracker_list) {
    std::list<std::list<PileupTracker*>> groups;
    GroupByAlignmentStart(*tracker_list, &groups);
    int num_discarded_items = LevelGroups(&groups);
    tracker_list->clear();
    for (const auto& group : groups) {
      tracker_list->insert(tracker_list->end(), group.begin(), group.end());
    }

    return num_discarded_items;
  }

  int GetToCoverage() const {
    return target_size_;
  }

 private:

  /**
   * Group the underlying readStatesByAlignmentStart into a list of list of alignment state machines,
   * where each list contains machines with a unique genome site.  The outer list is ordered
   * by alignment start.
   *
   * For example, if the flat list has alignment starts [10, 10, 11, 12, 12, 13] then
   * the resulting grouping will be [[10, 10], [11], [12, 12], [13]].
   *
   */
  void GroupByAlignmentStart(
      const std::list<PileupTracker*>& tracker_list,
      std::list<std::list<PileupTracker*>>* group_list) {
    PileupTracker* last = nullptr;
    for (const auto& item : tracker_list) {
      if (last == nullptr ||
          item->GetGenomeOffset() != last->GetGenomeOffset()) {
        // we've advanced to a place where the state machine has a different state,
        // so start a new list
        group_list->push_back({});
        last = item;
      }
      group_list->back().push_back(item);
    }
  }

  int LevelGroups(std::list<std::list<PileupTracker*>>* groups) {
    std::vector<int> group_sizes;
    group_sizes.reserve(groups->size());

    int total_size = 0;
    for (const auto& item : *groups) {
      group_sizes.push_back(item.size());
      total_size += item.size();
    }
    if (total_size < target_size_) {
      return 0;
    }

    // We will try to remove exactly this many items, however we will refuse to
    // allow any one group to fall below size 1, and so might end up
    // removing fewer items than this
    int num_items_to_remove = total_size - target_size_;

    int current_group_idx = 0;
    int num_consecutive_unmodifiable_groups = 0;

    // continue until we've either remove all the items we wanted to, or we
    // can't remove any more items without violating the constraint that
    // all groups must be left with at least one item
    while (num_items_to_remove > 0 &&
           num_consecutive_unmodifiable_groups < group_sizes.size()) {
      if (group_sizes[current_group_idx] > min_element_per_stack_) {
        group_sizes[current_group_idx]--;
        num_items_to_remove --;
        num_consecutive_unmodifiable_groups = 0;
      } else {
        num_consecutive_unmodifiable_groups ++;
      }

      current_group_idx = (current_group_idx + 1) % group_sizes.size();
    }

    // Now we actually go through and reduce each group to its new count as
    // specified in group_sizes
    current_group_idx = 0;
    int num_discarded_items = 0;
    for (auto iter = groups->begin(); iter != groups->end(); ++iter) {
      num_discarded_items += DownsampleOneGroup(
          group_sizes[current_group_idx++], &(*iter));
    }
    return num_discarded_items;
  }

  int DownsampleOneGroup(int num_items_to_keep, std::list<PileupTracker*>* group) {
    int group_size = group->size();
    if (num_items_to_keep >= group_size) return 0;

    std::vector<bool> item_to_keep(group_size, false);
    item_to_keep.reserve(group_size);
    std::vector<int> chosen_index = utils::SampleIndicesWithoutReplacemement(
        group_size,
        num_items_to_keep);
    for (int selected_index : chosen_index) {
      item_to_keep[selected_index] = true;
    }

    int num_discarded_items = 0;
    int current_index = 0;
    for (auto iter = group->begin(); iter != group->end();) {
      if (!item_to_keep[current_index++]) {
        iter = group->erase(iter);
        num_discarded_items ++;
      } else {
        ++iter;
      }
    }

    return num_discarded_items;
  }

  int target_size_;
  int min_element_per_stack_;

};

} // easehts
} // ncic

#endif
