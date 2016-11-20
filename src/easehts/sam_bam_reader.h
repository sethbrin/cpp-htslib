//
// Created by zp on 8/24/16.
//

#ifndef EASEHTSLIB_READER_H
#define EASEHTSLIB_READER_H

#include "utils.h"
#include "noncopyable.h"
#include "sam_bam_record.h"
#include "sam_sequence_dictionary.h"

#include <string>
#include <queue>

namespace ncic {

namespace easehts {

/**
 * The base class of sam/bam reader file
 * @example
 *     SAMBAMNormalReader reader(fp);
 *     SAMBAMRecord record;
 *     while (reader.HasNext(&record)) {
 *         // process the record
 *     }
 */
class AbstractSAMBAMReader : public NonCopyable {
 public:
  explicit AbstractSAMBAMReader(samFile* fp)
    : fp_(fp), header_(fp) {}

  explicit AbstractSAMBAMReader(const std::string& filename) {
    ERROR_COND(!utils::FileExists(filename),
        utils::StringFormatCStr("%s file is not exist!", filename.c_str()));
    samFile* fp = sam_open(filename.c_str(), "r");
    fp_ = fp;
    header_ = SAMSequenceDictionary(fp);
  }

  virtual bool HasNext(SAMBAMRecord* record) = 0;

  // load the bam data to record
  // just for the situation that use the raw data
  virtual bool HasNext(bam1_t* record) = 0;

  const SAMSequenceDictionary& GetSequenceDictionary() const {
    return header_;
  }

 protected:
  samFile* fp_;
  SAMSequenceDictionary header_;
};

typedef struct FileType {
  std::string name; // typename such as SAM
  std::string file_extension; // file extension .sam
  std::string index_extension; // index extension .bai

  const static FileType SRA_TYPE;
  const static FileType CRAM_TYPE;
  const static FileType BAM_TYPE;
  const static FileType SAM_TYPE;

} FileType;

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
class SAMBAMTextReader : public AbstractSAMBAMReader {
 public:
  using AbstractSAMBAMReader::AbstractSAMBAMReader;
  bool HasNext(SAMBAMRecord* record) override;
  bool HasNext(bam1_t* record) override;
};


class BAMIndex : public NonCopyable {
 public:
  explicit BAMIndex(const std::string& filename)
      : filename_(filename) {
    hts_idx_ = ::hts_idx_load(filename_.c_str(), HTS_FMT_BAI);
  }

  BAMIndex(const std::string& filename, const std::string& ext)
      : filename_(filename),
        filename_ext_(ext) {
    hts_idx_ = ::hts_idx_load2(filename_.c_str(), filename_ext_.c_str());
  }

  std::string GetFileName() {
    return filename_;
  }

  std::string GetFileExtName() {
    return filename_ext_;
  }

  hts_idx_t* GetRawIndex() {
    return hts_idx_;
  }

  ~BAMIndex() {
    ::hts_idx_destroy(hts_idx_);
  }

 private:
  std::string filename_; // the data filename
  std::string filename_ext_; // the default is .bai/.csi
  hts_idx_t* hts_idx_;
};

class BAMIndexReader : public AbstractSAMBAMReader {
 public:
  explicit BAMIndexReader(const std::string& filename)
      : AbstractSAMBAMReader(filename),
      index_(filename),
      hts_iter_(nullptr) {}

  bool HasNext(SAMBAMRecord* record) override;
  bool HasNext(bam1_t* record) override;

  /**
   * Set the region that need to read
   *
   * XXX Set a new region will clear the last region hts_itr info,
   * so you should read all the info, then set another region to read.
   *
   * @param region
   */
  void SetRegion(const std::string& region);
  void SetRegion(int32_t tid, int32_t begin, int32_t end);

 private:
  // bam region iters
  hts_itr_t* hts_iter_;
  BAMIndex index_;

};

/**
 * read many region in once
 */
class BAMIndexBatchReader : public AbstractSAMBAMReader {
 public:
  explicit BAMIndexBatchReader(const std::string& filename)
      : AbstractSAMBAMReader(filename),
      index_(filename) {}

  bool HasNext(SAMBAMRecord* record) override;
  bool HasNext(bam1_t* record) override;

  void AddRegion(const std::string& region);
  void AddRegion(int32_t tid, int32_t begin, int32_t end);

 private:
  // bam region iters
  std::queue<hts_itr_t*> hts_iters_;
  BAMIndex index_;

};




} // easehts
} // ncic

#endif //EASEHTSLIB_READER_H
