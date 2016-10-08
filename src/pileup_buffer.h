//
// Created by zp on 9/13/16.
//

#ifndef EASEHTSLIB_PILEUP_BUFFER_H
#define EASEHTSLIB_PILEUP_BUFFER_H

#include <list>
#include <mutex>
#include <atomic>
#include <memory>
#include <condition_variable>
#include <vector>
#include <htslib/vcf.h>
#include <map>

namespace ncic {
namespace easehts {
namespace mpileup {

typedef struct PackData {
  int idx;
  int tid;
  int begin;
  int end;
  volatile bool is_done;
  std::vector<bcf1_t*> data;
  PackData(int idx_, int tid_, int beg_, int end_)
      : idx(idx_),
        tid(tid_),
        begin(beg_),
        end(end_),
        is_done(false)
  { }
} PackData;

// pileup write buffer
class PileupBuffer {
 public:
  explicit PileupBuffer() : is_stop_(false), is_write_stop_(false), current_block_idx(0) {
  }

  void clearBuffer(PackData* pd) {
    for (bcf1_t* item : pd->data) {
      bcf_destroy1(item);
    }
    delete pd;
  }

  void Stop() {
    std::unique_lock<std::mutex> lock(empty_mutex_);
    is_stop_ = true;
    empty_cond_.notify_all();
  }

  void Add(int idx, int tid, int beg, int end) {
    assert(is_stop_ != true);
    PackData* pd = new PackData(idx, tid, beg, end);
    std::unique_lock<std::mutex> lock(empty_mutex_);
    assert(is_stop_ != true);
    tasks_.push_back(pd);
    data_map_[idx] = std::prev(tasks_.end());
    empty_cond_.notify_one();

    fprintf(stderr, "add region %d: %d:%d-%d\n", idx, tid, beg, end);
  }

  // Get the needed processed data iter
  PackData* GetProcessPackData(int idx) {
    std::unique_lock<std::mutex> lock(empty_mutex_);

    if (tasks_.empty() && is_stop_) {
      return NULL;
    }

    while (tasks_.empty()) {
      empty_cond_.wait(lock);
    }

    assert(data_map_.find(idx) != data_map_.end());
    return *(data_map_[idx]);
  }

  void SetPackDataDone(int idx) {
    assert(data_map_.find(idx) != data_map_.end());

    auto pd_iter = data_map_[idx];
    (*pd_iter)->is_done = true;

    fprintf(stderr, "region done %d: %d:%d-%d\n", idx, (*pd_iter)->tid, (*pd_iter)->begin, (*pd_iter)->end);
    std::unique_lock<std::mutex> lock(empty_mutex_);
    std::unique_lock<std::mutex> lock2(write_mutex_);

    if (current_block_idx < idx) {
      return;
    }

    while (pd_iter != tasks_.end()) {
      if ((*pd_iter)->is_done) {
        write_tasks_.push_back(*pd_iter);
        pd_iter = tasks_.erase(pd_iter);
        current_block_idx++;
      } else {
        break;
      }
    }
    // when tasks queue empty and the stop add packdata means the write buffer stop
    if (tasks_.empty() && is_stop_) {
      is_write_stop_ = true;
    }
    write_cond_.notify_one();
    fprintf(stderr, "notify: %d\n", is_stop_);
  }

  PackData* PopFront() {
    std::unique_lock<std::mutex> lock(write_mutex_);

    if (is_write_stop_ && write_tasks_.empty()) {
      return NULL;
    }

    while (write_tasks_.empty()) {
      fprintf(stderr, "wait: %s\n", is_write_stop_ ? "true":"false");
      write_cond_.wait(lock);
    }

    PackData* res = write_tasks_.front();
    fprintf(stderr, "write region done %d: %d:%d-%d\n", res->idx, res->tid, res->begin, res->end);
    data_map_.erase(res->idx);
    write_tasks_.pop_front();
    return res;
  }

  ~PileupBuffer() {
    is_stop_ = true;

    empty_cond_.notify_all();
  }
 private:
  std::list<PackData*> tasks_;
  std::list<PackData*> write_tasks_;
  int current_block_idx;
  std::map<int, std::list<PackData*>::iterator> data_map_;
  std::mutex empty_mutex_;
  std::mutex write_mutex_;
  std::condition_variable empty_cond_; // the tasks is not empty
  std::condition_variable write_cond_;
  bool is_stop_;
  bool is_write_stop_;
};

} // mpileup
} // easehts
} // ncic


#endif //EASEHTSLIB_PILEUP_BUFFER_H
