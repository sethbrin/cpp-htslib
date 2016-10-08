//
// Created by zp on 9/3/16.
//

#ifndef EASEHTSLIB_BAM_READER_H
#define EASEHTSLIB_BAM_READER_H

#include "bam_index.h"
#include "reader.h"

#include <memory>
#include <queue>

namespace ncic {

namespace easehts {

class BAMIndexReader : public AbstractReader {
 public:
  explicit BAMIndexReader(samFile* fp, const std::shared_ptr<BAMIndex>& index)
      : AbstractReader(fp),
        index_(index){}

  bool HasNext(SAMBAMRecord* record);

  void AddRegion(const char* region);
  void AddRegion(int32_t tid, int32_t begin, int32_t end);

 private:
  // bam region iters
  std::queue<std::shared_ptr<hts_itr_t>> hts_iters_;
  std::shared_ptr<BAMIndex> index_;

};

} // easehts
} // ncic


#endif //EASEHTSLIB_BAM_READER_H
