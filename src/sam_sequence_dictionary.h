//
// Created by zp on 8/24/16.
//

#ifndef EASEHTSLIB_SAM_SEQUENCE_DICTIONARY_H
#define EASEHTSLIB_SAM_SEQUENCE_DICTIONARY_H

#include "sam_seqence_record.h"

#include <vector>
#include <map>
#include <memory>

#include <assert.h>

namespace ncic {

namespace easehts {

/**
 * Collection of SAMSequenceRecords.
 */
class SAMSequenceDictionary {

 public:
  void AddSequenceRecord(const SAMSequenceRecord& sequence_record) {
    sequences_.push_back(sequence_record);
    sequences_map_[sequence_record.getSeqenceName()] = static_cast<int32_t>(sequences_.size()) - 1;
  }

  // XXX you should first decide if hasSequence
  SAMSequenceRecord& GetSequence(const std::string& name) {
    assert(HasSequence(name));

    return sequences_[sequences_map_[name]];
  }

  bool HasSequence(const std::string& name) const {
    if (sequences_map_.find(name) != sequences_map_.end()) {
      return true;
    }

    return false;
  }

  size_t Size() const {
    return sequences_.size();
  }

  SAMSequenceRecord& operator[](size_t idx) {
    return sequences_[idx];
  }

  SAMSequenceRecord& At(size_t idx) {
    assert(idx < sequences_.size());

    return sequences_[idx];
  }

 private:
  std::vector<SAMSequenceRecord> sequences_;
  // sequence name -> SAMSequenceRecord index
  std::map<std::string, int32_t> sequences_map_;
};

typedef std::shared_ptr<SAMSequenceDictionary> SAMSequenceDictionaryPtr;

} // easehts
} // ncic

#endif //EASEHTSLIB_SAM_SEQUENCE_DICTIONARY_H
