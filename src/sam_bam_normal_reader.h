//
// Created by zp on 8/24/16.
//

#ifndef EASEHTSLIB_SAMREADER_H
#define EASEHTSLIB_SAMREADER_H

#include "reader.h"
#include "sam_sequence_dictionary.h"
#include "sam_bam_record.h"
#include "sam_bam_header.h"
#include "noncopyable.h"

#include <string>
#include <map>
#include <htslib/sam.h>


namespace ncic {

namespace easehts {

/**
 * A normal reader for sam/bam file, no need bam index file or region
 * read all the data in the file
 *
 * @example
 *     SAMBAMNormalReader reader(fp);
 *     SAMBAMRecord record;
 *     while (reader.HasNext(&record)) {
 *         // process the record
 *     }
 */
class SAMBAMNormalReader : public AbstractReader, public NonCopyable {
 public:
  // Prepare to read a SAM text file.
  explicit SAMBAMNormalReader(samFile* fp) : AbstractReader(fp) {
    CHECK_NOTNULL(fp);

    ReadHeader();
  }

  const SAMBAMHeader& GetFileHeader() const {
    return AbstractReader::header_;
  }

  bool HasNext(SAMBAMRecord* record);

 private:
  void ReadHeader();

  /**
   * parse the sam record line and store to the SAMBAMRecord
   *
   * @param
   *     record: SAMBAMRecord
   */
  //void ParseRecordLine(SAMBAMRecord* record);

 private:
  //SAMFileHeaderPtr file_header_;

};


} // easehts
} // ncic

#endif //EASEHTSLIB_SAMREADER_H
