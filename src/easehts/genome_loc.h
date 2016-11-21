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
  GenomeLoc(const std::string& contig, int contig_index,
            int start, int stop)
    : contig_name_(contig),
    contig_index_(contig_index),
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
   * compare this genomeLoc's contig to another genome loc
   * @param that the genome loc to compare contigs with
   * @return 0 if equal, -1 if that.contig is greater, 1 if this contig is greater
   */
  int CompareContigs(const GenomeLoc& that) const {
    if (contig_index_ == that.contig_index_) {
      return 0;
    } else if (contig_index_ > that.contig_index_) {
      return 1;
    } else {
      return -1;
    }
  }

  std::string GetContig() const {
    return contig_name_;
  }

  int GetContigIndex() const {
    return contig_index_;
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
    return GenomeLoc(GetContig(), GetContigIndex(), GetStart(), GetStart());
  }

  GenomeLoc GetStopLocation() {
    return GenomeLoc(GetContig(), GetContigIndex(), GetStop(), GetStop());
  }

  std::string ToString() const {
    if (GenomeLoc::IsUnmapped(*this)) return "unmapped";
    if (ThroughEndOfContigP() && AtBeginningContigP()) {
      return GetContig();
    } else if (ThroughEndOfContigP() || GetStart() == GetStop()) {
      return utils::StringFormat("%s:%d", GetContig().c_str(), GetStart());
    } else {
      return utils::StringFormat("%s:%d-%d", GetContig().c_str(), GetStart(), GetStop());
    }

  }

  int contig_index_;
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
  int GetContigIndex(const std::string& contig) const;

 private:
    const SAMSequenceDictionary& ref_header_;
    ValidationLevel validation_level_;
};

class IntervalUtils : public NonCopyable {
 public:
  static std::vector<GenomeLoc> IntervalFileToList(const GenomeLocParser& gl_parser,
                                                   const std::string& filename);
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
