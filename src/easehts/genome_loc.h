//
// Created by zp on 9/3/16.
//

#ifndef EASEHTSLIB_GENOMELOC_H
#define EASEHTSLIB_GENOMELOC_H

#include "noncopyable.h"
#include "utils.h"
#include "sam_sequence_dictionary.h"

#include <limits.h>
#include <string>
#include <vector>

namespace ncic {

namespace easehts {

/**
 * Genome location representation.  It is *** 1 *** based closed.  Note that
 * GenomeLocs start and stop values can be any positive or negative number,
 * by design.  Bound validation is a feature of the GenomeLocParser,
 * and not a fundamental constraint of the GenomeLoc
 */
class GenomeLoc {
 public:
  /**
   * the basic components of a genome loc, its contig index,
   * start and stop position, and (optionally) the contig name
   */
  GenomeLoc(const std::string& contig, int contig_id,
            int start, int stop)
    : contig_name_(contig),
    contig_id_(contig_id),
    start_(start),
    stop_(stop) {}

  // WARN As contig_name_ may be expensive, here we just use contig_id
  // to construct it, but we can not get the ContigName
  // just for the temporary use
  GenomeLoc(int contig_id,
            int start, int stop)
    : contig_id_(contig_id),
    start_(start),
    stop_(stop) {}

  GenomeLoc(const std::string& contig)
    : GenomeLoc(contig, -1, 0, 0) {}

  GenomeLoc()
    : GenomeLoc(std::string(), -1, 0, 0) {}

  bool operator < (const GenomeLoc& that) const {
    return CompareTo(that) == -1;
  }

  bool operator > (const GenomeLoc& that) const {
    return CompareTo(that) == 1;
  }

  bool operator== (const GenomeLoc& that) const {
    return CompareTo(that) == 0;
  }

  /**
   * Test whether this contig is completely after contig 'that'
   * @param that contig to test against
   * @return true if this contig starts after 'that' end, false otherwise
   */
  bool IsPast(const GenomeLoc& that) const {
    int comparsion = CompareContigs(that);
    return comparsion == 1 ||
      (comparsion == 0 && start_ > that.stop_);
  }

  /**
   * Test whether this contig is completely before contig 'that'
   * @param that contig to test against
   * @return true if this contig starts before 'that' end, false otherwise
   */
  bool IsBefore(const GenomeLoc& that) const {
    int comparsion = CompareContigs(that);
    return comparsion == -1 ||
      (comparsion == 0 && stop_ < that.start_);
  }


  /**
   * compare this genomeLoc's contig to another genome loc
   * @param that the genome loc to compare contigs with
   * @return 0 if equal, -1 if that.contig is greater, 1 if this contig is greater
   */
  int CompareContigs(const GenomeLoc& that) const {
    if (contig_id_ == that.contig_id_) {
      return 0;
    } else if (contig_id_ > that.contig_id_) {
      return 1;
    } else {
      return -1;
    }
  }

  std::string GetContig() const {
    return contig_name_;
  }

  int GetContigId() const {
    return contig_id_;
  }

  int GetStart() const {
    return start_;
  }

  int GetStop() const {
    return stop_;
  }

  GenomeLoc GetLocation() {
    return *this;
  }

  GenomeLoc GetStartLocation() {
    return GenomeLoc(GetContig(), GetContigId(), GetStart(), GetStart());
  }

  GenomeLoc GetStopLocation() {
    return GenomeLoc(GetContig(), GetContigId(), GetStop(), GetStop());
  }

  bool DisjointP(const GenomeLoc& that) const {
    return contig_id_ != that.contig_id_ ||
      start_ > that.stop_ ||
      that.start_ > stop_;
  }

  bool DisContinuousP(const GenomeLoc& that) const {
    return contig_id_ != that.contig_id_ ||
      (start_ - 1) > that.stop_ ||
      (that.start_ - 1) > stop_;
  }

  bool OverlapsP(const GenomeLoc& that) const {
    return !DisjointP(that);
  }

  bool ContiguousP(const GenomeLoc& that) const {
    return !DisContinuousP(that);
  }

  GenomeLoc Merge(const GenomeLoc& that) const {
    if (GenomeLoc::IsUnmapped(*this) ||
        GenomeLoc::IsUnmapped(that)) {
      ERROR_COND(!GenomeLoc::IsUnmapped(*this) ||
                 !GenomeLoc::IsUnmapped(that),
                 utils::StringFormatCStr(
                     "Trid to merge a mapped and an unmapped genome loc"));
      return kUnmapped;
    }

    ERROR_COND(!ContiguousP(that),
               utils::StringFormatCStr(
                   "The two genome loc's need to be contiguous"));

    return GenomeLoc(GetContig(), contig_id_,
                     std::min(GetStart(), that.GetStart()),
                     std::max(GetStop(), that.GetStop()));
  }

  std::string ToString() const {
    if (GenomeLoc::IsUnmapped(*this)) return "unmapped";
    if (ThroughEndOfContigP() && AtBeginningContigP()) {
      return GetContig();
    } else if (ThroughEndOfContigP() || GetStart() == GetStop()) {
      return utils::StringFormat("%s:%d", GetContig().c_str(), GetStart());
    } else {
      return utils::StringFormat("%s:%d-%d",
                                 GetContig().c_str(), GetStart(), GetStop());
    }

  }

  int contig_id_;
  int start_;
  int stop_;
  std::string contig_name_;

  static const GenomeLoc kUnmapped;
  static const GenomeLoc kWholeGenome;

  static bool IsUnmapped(const GenomeLoc& loc) {
    return &loc == &kUnmapped;
  }


 private:
  bool ThroughEndOfContigP() const { return stop_ == INT_MAX; }
  bool AtBeginningContigP() const { return start_ == 1; }

  int CompareTo(const GenomeLoc& that) const {
    int result = 0;
    if (this == &that) {
      result = 0;
    } else if (GenomeLoc::IsUnmapped(*this)) {
      result = 1;
    } else if (GenomeLoc::IsUnmapped(that)) {
      result = -1;
    } else {
      int cmp_contig = CompareContigs(that);

      if (cmp_contig != 0) {
        result = cmp_contig;
      } else {
        if (GetStart() < that.GetStart()) {
          result = -1;
        } else if (GetStart() > that.GetStart()) {
          result = 1;
        } else if (GetStop() < that.GetStop()) {
          result = -1;
        } else if (GetStop() > that.GetStop()) {
          result = 1;
        }
      }
    }
    return result;
  }

};



/**
 * Factory class for creating GenomeLocs
 */
class GenomeLocParser : public NonCopyable {
  /*
   * How much validation should we do at runtime with this parser?
   */
  enum ValidationLevel {
    /** Do the standard amount of validation */
    STANDARD,
    /** Don't do any real checking at all */
    NONE
  };

 public:
  GenomeLocParser(const SAMSequenceDictionary& header, ValidationLevel level)
    : ref_header_(header),
    validation_level_{level} {}

  explicit GenomeLocParser(const SAMSequenceDictionary& header)
    : ref_header_(header),
    validation_level_{ValidationLevel::STANDARD} {}

  /**
   * Would a genome loc created with a given parameters be valid w.r.t. the
   * master sequence dictionary?
   * @param contig the contig we'd use
   * @param start the start position
   * @param stop the stop
   * @param must_be_on_reference should we require the resulting genome loc to
   * be completely on rhe reference
   * @return true if this would produce a valid genome loc, false otherwise
   */
  bool ValidateGenomeLoc(const std::string& contig, int start,
                         int stop, bool must_be_on_reference) const;
  int GetContigId(const std::string& contig) const;

 private:
    const SAMSequenceDictionary& ref_header_;
    ValidationLevel validation_level_;
};

class IntervalUtils : public NonCopyable {
 public:
  static std::vector<GenomeLoc> LoadIntervals(
      const GenomeLocParser& gl_parser,
      const std::string& filename,
      const int interval_padding);

  static std::vector<GenomeLoc> IntervalFileToList(
      const GenomeLocParser& gl_parser,
      const std::string& filename,
      const int interval_padding);

  static std::vector<GenomeLoc> MergeIntervals(
      const std::vector<GenomeLoc>& genome_locs);

  /**
   * read interval header from filename
   * @param filename the filename contains header, such as sam/interval
   * @return SAMSequenceDictionary
   */
  static SAMSequenceDictionary ReadHeader(htsFile* fp);

};

} // easehts
} // ncic

#endif
