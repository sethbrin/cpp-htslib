//
// Created by zp on 12/15/16.
//

#ifndef EASEHTSLIB_VCF_H_
#define EASEHTSLIB_VCF_H_

#include "easehts/utils.h"
#include "easehts/genome_loc.h"

#include <htslib/sam.h>
#include <htslib/vcf.h>
#include <htslib/tbx.h>

#include <cassert>
#include <list>
#include <memory>

namespace ncic {
namespace easehts {

class VCFHeader {
 public:
  explicit VCFHeader(const char *mode) {
    header_ = ::bcf_hdr_init(mode);
  }

  explicit VCFHeader(bcf_hdr_t* h)
    : header_(h) {}

  ~VCFHeader() {
    ::bcf_hdr_destroy(header_);
  }

  bcf_hdr_t* GetRawHeader() const {
    return header_;
  }

  std::string GetVersion() const {
    return ::bcf_hdr_get_version(header_);
  }

  int GetSamplesCnt() const {
    return bcf_hdr_nsamples(header_);
  }

  int GetContigId(const std::string& contig) const {
    return ::bcf_hdr_name2id(header_, contig.c_str());
  }

  std::string GetContig(int id) const {
    return ::bcf_hdr_id2name(header_, id);
  }

  int GetContigLength(const std::string& contig) const {
    return header_->id[BCF_DT_CTG][GetContigId(contig)].val->info[0];
  }

  int GetContigLength(int id) const {
    return header_->id[BCF_DT_CTG][id].val->info[0];
  }

  GenomeLoc CreateOverEntireContig(int id) const {
    return GenomeLoc(id, 0, GetContigLength(id) - 1);
  }

 private:
  bcf_hdr_t* header_;
};


class VariantContext {
 public:
  VariantContext() {
    record_ = ::bcf_init1();
  }

  ~VariantContext() {
    ::bcf_destroy1(record_);
  }

  bcf1_t* GetRawRecord() const {
    return record_;
  }

  int GetId() const {
    return record_->rid;
  }

  int GetPos() const {
    return record_->pos;
  }

  int GetReferenceLength() const {
    return record_->rlen;
  }

  float GetQual() const {
    return record_->qual;
  }

  const GenomeLoc& GetLocation() const {
    if (location_ == nullptr) {
      location_.reset(new GenomeLoc(GetId(), GetPos(), GetPos()));
    }
    return *location_;
  }

 private:
  bcf1_t* record_;

  mutable std::unique_ptr<GenomeLoc> location_;
};

class VCFReader {
 public:
  explicit VCFReader(const std::string& filename,
                     const std::string& mode="r") {
    ERROR_COND(!utils::FileExists(filename),
        utils::StringFormatCStr("%s file is not exist!", filename.c_str()));
    filename_ = filename;
    fp_ = ::hts_open(filename.c_str(), mode.c_str());
    header_.reset(new VCFHeader(::bcf_hdr_read(fp_)));
    cur_record_ = nullptr;
  }

  ~VCFReader() {
    int ret;
    if ((ret = ::hts_close(fp_))) {
      fprintf(stderr, "hts_close(%s) non zero status %d\n",
              filename_.c_str(), ret);
      ::exit(ret);
    }
  }

  virtual bool HasNext() = 0;

  VariantContext* Next() {
    return cur_record_;
  }

  const VCFHeader& GetHeader() const {
    return *header_;
  }

 protected:
  std::unique_ptr<VCFHeader> header_;
  std::string filename_;
  htsFile* fp_;

  VariantContext* cur_record_;
};

class VCFTextReader : public VCFReader {
 public:
  using VCFReader::VCFReader;

  bool HasNext() override {
    VariantContext* record = new VariantContext();

    if (::bcf_read1(fp_, header_->GetRawHeader(),
                    record->GetRawRecord()) >= 0) {
      cur_record_ = record;
      return true;
    }
    delete cur_record_;
    return false;
  }

};

class VCFIndexReader : public VCFReader {
 public:
  explicit VCFIndexReader(const std::string& filename,
                          const std::string& mode="r")
    : VCFReader(filename, mode) {
    tbx_idx_ = tbx_index_load(filename_.c_str());
    WARN_COND(tbx_idx_ == nullptr,
              utils::StringFormatCStr(
                  "the index of file %s can not loaded", filename.c_str()));
    hts_iter_ = nullptr;
  }

  ~VCFIndexReader() {
    if (hts_iter_) tbx_itr_destroy(hts_iter_);
    if (tbx_idx_) tbx_destroy(tbx_idx_);
  }

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


  bool HasNext() override;

 private:
  tbx_t* tbx_idx_;
  hts_itr_t* hts_iter_;
};

class VCFTraverse {
 public:
  VCFTraverse(VCFReader* reader)
    : reader_(reader)
  {}

  ~VCFTraverse() {
    for (auto iter = buffer_list_.begin();
         iter != buffer_list_.end(); ++iter) {
      delete *iter;
    }
    buffer_list_.clear();
  }

  VCFReader* GetReader() const {
    return reader_;
  }

  /**
   * Seeks forward through the file until the specified interval is reached.
   * The location object <code>interval</code> can be either a single point or
   * an extended interval. All records overlapping with the whole interval
   * will be returned, or null if no such records exist.
   * Query interval must start at or after the iterator's current location,
   * or exception will be thrown.
   *
   * Query interval must end at or after the stop position of the previous
   * query, if any, or an error be triggered: subsequent queries that
   * end before the stop of previous ones are illegal.
   *
   * If SeekFroward() is called, we can call GetFirstRecord()
   * OR GetRecords() to get the VariantContext associated with the interval
   *
   * @example
   * SeekFroward(interval);
   * VariantContext* record = GetFirstRecord();
   *
   * @param interval point-like genomic location to fastforward to.
   * @return ROD object at (or overlapping with) the specified position, or null if no such ROD exists.
   */
  void SeekFroward(const GenomeLoc& interval);

  VariantContext* GetFirstRecord();

  void GetRecords(std::vector<VariantContext*>* records);
 private:

  /**
   * Removes records that end before the curr_position from the list of
   * currently kept records. This is a
   * convenience (private) shortcut that does not perform extensive
   * checking. In particular, it assumes that
   * curr_position <= max_position, as well as that we are
   * still on the same contig.
   */
  void PurgeOutofScopeRecords();

  VCFReader* reader_;

  // here we will keep a pile of records overlaping with current position; when
  // we traverse it
  std::list<VariantContext*> buffer_list_;

  // current traverse record
  VariantContext* cur_record_ = nullptr;

  int cur_contig_id_ = -1;
  // where the iterator is currently positioned on the genome
  int cur_pos_ = 0;
  // the rightmost stop position of currently loaded records
  int max_pos_ = 0;

  // the stop position of the last query. We can query only in forward
  // direction ("seek forward"); it is not only the start position of every
  // successive query that can not be before the start of the previous one
  // (curr_start), but it is also illegal for a query interval to *end* before
  // the end of previous query, otherwise we can end up in an inconsistent
  // state
  int cur_query_end_ = -1;
};

} // easehts
} // ncic

#endif
