//
// Created by zp on 8/24/16.
//

#ifndef EASEHTSLIB_SAM_SEQENCE_RECORD_H
#define EASEHTSLIB_SAM_SEQENCE_RECORD_H


#include <string>

namespace ncic {

namespace easehts {

/*
 * Header information about a reference sequence.  Corresponds to @SQ header record in SAM text header.
 */
class SAMSequenceRecord {
 public:
  explicit SAMSequenceRecord(const std::string& name)
      : sequence_name_(name),
        sequence_length_(UNKNOWN_SEQUENCE_LENGTH) {
  }

  SAMSequenceRecord(const std::string& name, int len)
      : sequence_name_(name),
        sequence_length_(len) {
  }

  const std::string& getSeqenceName() const {
    return sequence_name_;
  }

  int getSequenceLength() const {
    return sequence_length_;
  }

 private:
  std::string sequence_name_;
  int sequence_length_;
  const static int UNKNOWN_SEQUENCE_LENGTH = 0;
};

} // easehts
} // ncic

#endif //EASEHTSLIB_SAM_SEQENCE_RECORD_H
