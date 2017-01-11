//
// Created by zp on 12/15/16.
//

#ifndef EASEHTSLIB_VCF_H_
#define EASEHTSLIB_VCF_H_

#include "easehts/utils.h"
#include "easehts/genome_loc.h"
#include "easehts/noncopyable.h"
#include "easehts/sam_sequence_dictionary.h"
#include "easehts/genotype.h"

#include <htslib/sam.h>
#include <htslib/vcf.h>
#include <htslib/tbx.h>

#include <cassert>
#include <cmath>
#include <list>
#include <memory>
#include <mutex>
#include <queue>
#include <set>
#include <string>
#include <unordered_map>

namespace ncic {
namespace easehts {

class VCFHeader : public NonCopyable {
 public:
  explicit VCFHeader(const char *mode) {
    header_ = ::bcf_hdr_init(mode);
  }
  explicit VCFHeader(const std::string& mode) {
    header_ = ::bcf_hdr_init(mode.c_str());
  }

  explicit VCFHeader(bcf_hdr_t* h)
    : header_(h) {}

  void Copy(const VCFHeader& rhs) {
    header_ = bcf_hdr_dup(rhs.header_);
  }

  ~VCFHeader() {
    ::bcf_hdr_destroy(header_);
  }

  bcf_hdr_t* GetRawHeader() const {
    return header_;
  }

  std::string GetVersion() const {
    return ::bcf_hdr_get_version(header_);
  }

  void SetVersion(const std::string& version) {
    bcf_hdr_set_version(header_, version.c_str());
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

  void AddHeaderLine(const std::string& line) {
    bcf_hdr_append(header_, line.c_str());
  }

  void AddContigs(bam_hdr_t* hdr, const std::string& assemby) {
    for (int idx = 0; idx < hdr->n_targets; idx++) {
      std::string line = "##contig=<";
      line += "ID=" + std::string(hdr->target_name[idx]) + ",";
      line += "length=" + std::to_string(hdr->target_len[idx]) + ",";
      // TODO assemby may change, current just set it to b37
      line += "assembly=" + assemby + ">";
      bcf_hdr_append(header_, line.c_str());
    }
  }

  void AddSample(const std::string& sample_name) {
    bcf_hdr_add_sample(header_, sample_name.c_str());
  }
  void FinishSample() {
    bcf_hdr_add_sample(header_, nullptr);
  }

 private:
  bcf_hdr_t* header_;
};


class VariantContext {
 public:
  VariantContext() {
    record_ = ::bcf_init1();
  }
  VariantContext(bcf1_t* record) {
    record_ = bcf_dup(record);
  }

  ~VariantContext() {
    ::bcf_destroy1(record_);
  }

  bcf1_t* GetRawRecord() const {
    return record_;
  }

  int GetContigId() const {
    return record_->rid;
  }

  void SetContigId(int id) {
    record_->rid = id;
  }

  void SetId(const std::string id, bcf_hdr_t* hdr) {
    bcf_update_id(hdr, record_, id.c_str());
  }

  std::string GetId() const {
    bcf_unpack(record_, BCF_UN_STR);
    return std::string(record_->d.id);
  }

  int GetPos() const {
    return record_->pos;
  }

  void SetPos(int pos) {
    record_->pos = pos;
  }

  int GetReferenceLength() const {
    return record_->rlen;
  }

  void SetReferenceLength(int len) {
    record_->rlen = len;
  }

  float GetQual() const {
    return record_->qual;
  }

  void Filter(const std::string& str, bcf_hdr_t* hdr) {
    int32_t tmpi = bcf_hdr_id2int(hdr, BCF_DT_ID, str.c_str());
    bcf_update_filter(hdr, record_, &tmpi, 1);
  }

  const GenomeLoc& GetLocation() const {
    if (location_ == nullptr) {
      location_.reset(new GenomeLoc(GetContigId(), GetPos(), GetPos()));
    }
    return *location_;
  }

  void SetLocation(const GenomeLoc& loc) {
    SetContigId(loc.GetContigId());
    SetPos(loc.GetStart());
  }

  void UpdateAlleleStr(const std::string& alleles, bcf_hdr_t* hdr) {
    bcf_update_alleles_str(hdr, record_, alleles.c_str());
  }

  template <typename T>
  void UpdateInfo(const std::string& key, T* val, int size, bcf_hdr_t* hdr) {
    ERROR("UnSupported value type");
  }

  void UpdateInfo(const std::string& key, int* val,
                       int size, bcf_hdr_t* hdr) {
    bcf_update_info_int32(hdr, record_, key.c_str(), val, size);
  }

  void UpdateInfo(const std::string& key, float* val,
                         int size, bcf_hdr_t* hdr) {
    bcf_update_info_float(hdr, record_, key.c_str(), val, size);
  }

  void UpdateStringInfo(const std::string& key,
                        const std::string& val, bcf_hdr_t* hdr) {
    bcf_update_info_string(hdr, record_, key.c_str(), val.c_str());
  }

  void UpdateFlagInfo(const std::string& key, bcf_hdr_t* hdr) {
    bcf_update_info_flag(hdr, record_, key.c_str(), nullptr, 1);
  }

  template <typename T>
  void UpdateFormat(const std::string& key, T* val,
                    int size, bcf_hdr_t* hdr) {
    ERROR("UnSupported value type");
  }

  void UpdateFormat(const std::string& key, int* val,
                    int size, bcf_hdr_t* hdr) {
    bcf_update_format_int32(hdr, record_, key.c_str(), val, size);
  }

  void UpdateFormat(const std::string& key, float* val,
                    int size, bcf_hdr_t* hdr) {
    bcf_update_format_float(hdr, record_, key.c_str(), val, size);
  }

  void UpdateGenetype(int* val, int size, bcf_hdr_t* hdr) {
    bcf_update_genotypes(hdr, record_, val, size);
  }

  static int GenotypePhaseInt(int num) {
    return bcf_gt_phased(num);
  }

  static int GenotypeUnPhaseInt(int num) {
    return bcf_gt_unphased(num);
  }

#define INT8_MISSING bcf_int8_missing
#define INT16_MISSING bcf_int16_missing
#define INT32_MISSING bcf_int32_missing
#define INT8_VECTOR_MISSING bcf_int8_vector_end
#define INT16_VECTOR_MISSING bcf_int16_vector_end
#define INT32_VECTOR_MISSING bcf_int32_vector_end
#define STR_MISSING bcf_str_missing
#define GT_MISSING bcf_gt_missing

 private:
  mutable bcf1_t* record_;

  mutable std::unique_ptr<GenomeLoc> location_;
  // TODO As current only support single char allele
  // They are all predefined in the Allele class,
  // so we can use address to check the equality
  // std::unordered_map<Allele*, int> allele_map;

};

class VCFConstants {
 public:
  const static std::string kAncestralAlleleKey;
  const static std::string kAlleleCountKey;
  const static std::string kMleAlleleCountKey;
  const static std::string kAlleleFrequenceKey;
  const static std::string kMleAlleleFrequenceKey;
  const static std::string kMlePerSampleAlleleCountKey;
  const static std::string kMlePerSampleAlleleFractionKey;
  const static std::string kAlleleNumberKey;
  const static std::string kRmsBaseQualityKey;
  const static std::string kCigarKey;
  const static std::string kDbsnpKey;
  const static std::string kDepthKey;
  const static std::string kDownsampleKey;
  const static std::string kExpectedAlleleCountKey;
  const static std::string kEndKey;
  const static std::string kGenotypeFilterKey;
  const static std::string kGenotypeKey;
  const static std::string kGenotypePosteriorsKey;
  const static std::string kGenotypeQualityKey;
  const static std::string kGenotypeAlleleDepths;
  const static std::string kGenotypePlKey;
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

  virtual ~VCFReader() {
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
    delete record;
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

  ~VCFIndexReader() override {
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
    if (cur_record_) delete cur_record_;
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
   * @return ROD object at (or overlapping with) the specified position,
   * or null if no such ROD exists.
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

// the class write variant context
class VariantContextWriter : public NonCopyable {
 public:
  VariantContextWriter(const std::string& filename,
                       const std::string& mode="w")
  : filename_(filename),
  header_("w") {
    fp_ = hts_open(filename.c_str(), mode.c_str());
    ERROR_COND(fp_ == nullptr, utils::StringFormatCStr(
            "Error open filename %s to write", filename.c_str()));
    is_closed_ = false;
  }

  ~VariantContextWriter() {
    Close();
  }

  VCFHeader* GetHeader() {
    return &header_;
  }

  virtual void WriteHeader() = 0;

  /**
   * Attempt to close vcf/bcf file
   */
  virtual void Close() {
    if (is_closed_) return;

    int ret = 0;
    if ((ret = ::hts_close(fp_))) {
      fprintf(stderr, "hts_close(%s) non zero status %d\n",
              filename_.c_str(), ret);
      ::exit(ret);
    }
    is_closed_ = true;
  }

  /**
   * Add a record to write
   */
  virtual void Add(std::unique_ptr<VariantContext>& vc) = 0;

 protected:
  htsFile* fp_;
  VCFHeader header_;
  bool is_closed_;
  std::string filename_;
};

class VCFWriter : public VariantContextWriter {
 public:
  using VariantContextWriter::VariantContextWriter;

  void WriteHeader() override {
    int ret = 0;
    {
      std::lock_guard<std::mutex> lock(mtx_);
      ret = bcf_hdr_write(fp_, header_.GetRawHeader());
    }
    ERROR_COND(ret < 0, "Failed to write header");
  }

  void Add(std::unique_ptr<VariantContext>& vc) override {
    int ret = 0;
    {
      std::lock_guard<std::mutex> lock(mtx_);
      ret = bcf_write(fp_, header_.GetRawHeader(), vc->GetRawRecord());
    }
    ERROR_COND(ret < 0, "Failed to write record");
  }

  void Add(std::unique_ptr<VariantContext>&& vc) {
    Add(vc);
  }

 private:
  std::mutex mtx_;
};

/**
 * This class writes VCF files, allowing records to be passed in unsorted.
 * It also enforces that it is never passed records of the same chromosome
 * with any other chromosome in between them.
 */
class SortingVariantContextWriter {
 public:
  SortingVariantContextWriter(VariantContextWriter* inner_writer,
                              int max_caching_start_distance,
                              bool take_ownership_of_inner)
    : inner_writer_(inner_writer),
    max_caching_start_distance_(max_caching_start_distance),
    take_ownership_of_inner_(take_ownership_of_inner) {
    most_upstream_writable_loc_ = kBeforeMostUpstreamLoc;
  }

  SortingVariantContextWriter(VariantContextWriter* inner_writer,
                              int max_caching_start_distance)
    : SortingVariantContextWriter(inner_writer,
                                  max_caching_start_distance, false) {}

  ~SortingVariantContextWriter() {
    Close();
  }

  void Close() {
    StopWaitingToSort();
    if (take_ownership_of_inner_) {
      inner_writer_->Close();
    }
  }

  VCFHeader* GetHeader() {
    return inner_writer_->GetHeader();
  }

  void NoteCurrentRecord(VariantContext* vc) {
    ERROR_COND(vc->GetPos() < most_upstream_writable_loc_,
               utils::StringFormatCStr(
                   "Permitted to write any record upstream of position %d, "
                   "but a record at %d:%d was just added",
                   most_upstream_writable_loc_, vc->GetContigId(),
                   vc->GetPos()));
    most_upstream_writable_loc_ = std::max(
        0,
        vc->GetPos() - max_caching_start_distance_);
  }

  void Add(std::unique_ptr<VariantContext>&& vc) {
    Add(vc);
  }
  /**
   * Add a record to the file
   *
   * @param vc the variant context
   */
  void Add (std::unique_ptr<VariantContext>& vc) {
    /**
     * NOTE that the code below does not prevent the successive
     * add()-ing of: (chr1, 10), (chr20, 200), (chr15, 100)
     * since there is no implicit ordering of chromosomes:
     */
    std::lock_guard<std::mutex> lock(mtx_);
    if (!queue_.empty()) {
      const std::unique_ptr<VariantContext>& first_record =
        queue_.top();
      if (vc->GetContigId() != first_record->GetContigId()) {
        ERROR_COND(finished_chromosomes_.find(vc->GetContigId()) !=
                   finished_chromosomes_.end(),
                   utils::StringFormatCStr(
                       "Add a record at %d:%d, but "
                       "already finished with chromosome",
                       vc->GetContigId(), vc->GetPos()));

        finished_chromosomes_.insert(first_record->GetContigId());
        StopWaitingToSort();
      }
    }

    NoteCurrentRecord(vc.get());
    queue_.push(std::move(vc));
    EmitSafeRecords();
  }

  void WriteHeader() {
    inner_writer_->WriteHeader();
  }

 private:
  void StopWaitingToSort() {
    EmitRecords(true);
    most_upstream_writable_loc_ = 0;
  }

  void EmitSafeRecords() {
    EmitRecords(false);
  }

  void EmitRecords(bool emit_unsafe) {
    while(!queue_.empty()) {
      const std::unique_ptr<VariantContext>& first_record =
        queue_.top();

      // No need to wait, waiting for nothing, or before what we're waiting for
      if (emit_unsafe || first_record->GetPos()
          <= most_upstream_writable_loc_) {
        inner_writer_->Add(
            const_cast<std::unique_ptr<VariantContext>&>(first_record));
        queue_.pop();
      } else {
        break;
      }
    }
  }

  class VariantContextCmp {
   public:
    bool operator()(const std::unique_ptr<VariantContext>& lhs,
                    const std::unique_ptr<VariantContext>& rhs) const {
      return lhs->GetPos() > rhs->GetPos();
    }
  };

  // The VCFWriter to which to actually write the sorted VCF records
  VariantContextWriter* inner_writer_;
  // the current queue of un-emitted records
  std::priority_queue<std::unique_ptr<VariantContext>,
    std::vector<std::unique_ptr<VariantContext>>, VariantContextCmp> queue_;

  // The locus until which we are permitted to write out(inclusive)
  int most_upstream_writable_loc_;
  // the maximum START distance between records that we'll cache
  int max_caching_start_distance_;
  const static int kBeforeMostUpstreamLoc;

  // The set of chromosomes already passed over and to which
  // it is forbidden to return
  std::set<int> finished_chromosomes_;

  // Should we call inner_writer_.close() in close
  bool take_ownership_of_inner_;

  std::mutex mtx_;
};

} // easehts
} // ncic

#endif
