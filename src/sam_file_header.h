//
// Created by zp on 8/24/16.
//

#ifndef EASEHTSLIB_SAM_FILE_HEADER_H
#define EASEHTSLIB_SAM_FILE_HEADER_H

#include "sam_sequence_dictionary.h"

#include <vector>
#include <memory>
#include <map>
#include <string>
#include <htslib/sam.h>

namespace ncic {

namespace easehts {

/**
 * Header information from a SAM or BAM file.
 */
class SAMFileHeader {
 public:
  SAMFileHeader() {
    sequence_dictionary_ = std::make_shared<SAMFileHeader>();
  }
  void SetSequenceDictionary(SAMFileHeaderPtr& sequence_dictionary) {
    sequence_dictionary_ = sequence_dictionary;
  }

  const SAMFileHeaderPtr& GetSequenceDictionary() const {
    return sequence_dictionary_;
  }

  void AddAttributes(const std::map<std::string, std::string> attributes) {
    attributes_.insert(attributes.begin(), attributes.end());
  }

 private:
  SAMFileHeaderPtr sequence_dictionary_;
  std::map<std::string, std::string> attributes_;

};

typedef std::shared_ptr<SAMFileHeader> SAMFileHeaderPtr;

} // easehts
} // ncic

#endif //EASEHTSLIB_SAM_FILE_HEADER_H
