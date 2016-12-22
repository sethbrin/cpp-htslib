//
// Created by zp on 12/21/16.
//

#ifndef EASEHTSLIB_GENOTYPE_H_
#define EASEHTSLIB_GENOTYPE_H_

#include <algorithm>
#include <string>
#include <memory>
#include <vector>
#include <unordered_map>

#include "easehts/utils.h"

namespace ncic {
namespace easehts {

class Allele {
 public:
  Allele(const std::string& bases, bool is_ref) {
    ERROR("Current UnSupported!");
    is_ref_ = false;
    is_no_call_ = false;
    is_symboic_ = false;

    ERROR_COND(WouldBeNullAllele(bases),
               "Null alleles are not supported");

    if (WouldBeNoCallAllele(bases)) {
      is_no_call_ = true;
      ERROR_COND(is_ref_, "Cannot tag a NoCall allele as the reference");
    } else {
      if (WouldBeSymbolicAllele(bases)) {
        is_symboic_ = true;
        ERROR_COND(is_ref_, "Cannot tag a symbolic allele as the reference");
      } else {
        bases_ = bases;
        std::transform(bases_.begin(), bases_.end(),
                       bases_.begin(), ::toupper);

        is_ref_ = is_ref;

        ERROR_COND(!AcceptableAlleleBases(bases),
                   "Unexpected base in allele bases");
      }
    }
  }

  /**
   * Create a new Allele that includes bases and if tagged as the reference
   * allele if isRef == true.  If bases
   * == '-', a Null allele is created.  If bases ==  '.',
   * a no call Allele is created. If bases ==  '*', a spanning deletions
   * Allele is created.
   *
   * @param bases the DNA sequence of this variation, '-', '.', or '*'
   * @param isRef should we make this a reference allele?
   */
  static Allele* Create(const std::string& bases, bool is_ref) {
    ERROR_COND(bases.empty(), "create: the Allele base string cannot be null; "
               "use new Allele() or new Allele(\"\") to create a Null allele");

    if ( bases.size() == 1 ) {
      return Create(bases[0], is_ref);
    } else {
      // TODO
      return new Allele(bases, is_ref);
    }
  }

  static Allele* Create(char base, bool is_ref) {
    // optimization to return a static constant Allele
    // for each single base object
    switch (base) {
      case '.':
        WARN_COND(is_ref,
                  "Cannot tag a NoCall allele as the reference allele");
        return kNoCall.get();
      case '*':
        WARN_COND(is_ref,
                  "Cannot tag a spanning deletions allele "
                  "as the reference allele");
        return kSpanDel.get();
      case 'A': case 'a' : return is_ref ? kRefA.get() : kAltA.get();
      case 'C': case 'c' : return is_ref ? kRefC.get() : kAltC.get();
      case 'G': case 'g' : return is_ref ? kRefG.get() : kAltG.get();
      case 'T': case 't' : return is_ref ? kRefT.get() : kAltT.get();
      case 'N': case 'n' : return is_ref ? kRefN.get() : kAltN.get();
      default:
        ERROR(utils::StringFormatCStr(
                "Illegal base [%c] seen in the allele", base));
    }
    return nullptr;
  }

  static Allele* Create(const std::string& bases) {
    return Create(bases, false);
  }

  static Allele* Create(char base) {
    return Create(base, false);
  }


  static bool AcceptableAlleleBases(const std::string& bases) {
    ERROR("UnImplemented function");
    return false;
  }

  static bool WouldBeNullAllele(const std::string& bases) {
    return (bases.size() == 1 && bases[0] == 45) ||
      bases.size() == 0;
  }

  static bool WouldBeNoCallAllele(const std::string& bases) {
    return bases.size() == 1 && bases[0] == 46;
  }

  static bool WouldBeSymbolicAllele(const std::string& bases) {
    if (bases.size() <= 2) return false;

    return bases[0] == 60 && bases.back() == 62 ||
      bases.find("[") != std::string::npos ||
      bases.find("]") != std::string::npos;
  }

 private:

  bool is_ref_;
  bool is_no_call_;
  bool is_symboic_;
  std::string bases_;

  const static std::string kNoCallString;
  const static std::string kSpanDelString;
  const static std::shared_ptr<Allele> kRefA;
  const static std::shared_ptr<Allele> kAltA;
  const static std::shared_ptr<Allele> kRefC;
  const static std::shared_ptr<Allele> kAltC;
  const static std::shared_ptr<Allele> kRefG;
  const static std::shared_ptr<Allele> kAltG;
  const static std::shared_ptr<Allele> kRefT;
  const static std::shared_ptr<Allele> kAltT;
  const static std::shared_ptr<Allele> kRefN;
  const static std::shared_ptr<Allele> kAltN;
  const static std::shared_ptr<Allele> kNoCall;
  const static std::shared_ptr<Allele> kSpanDel;
};

class GenotypeBuilder {
 public:
  GenotypeBuilder(const std::string& sample_name,
                  const std::vector<Allele*>& alleles)
    : sample_name_(sample_name),
    alleles_(alleles) {}

  GenotypeBuilder& DP(int DP) {
    this->DP_ = DP;
    return *this;
  }

  GenotypeBuilder& GQ(int GQ) {
    this->GQ_ = GQ;
    return *this;
  }

  GenotypeBuilder& AD(const std::vector<int>& AD) {
    this->AD_ = AD;
    return *this;
  }

  GenotypeBuilder& PL(const std::vector<int>& PL) {
    this->PL_ = PL;
    return *this;
  }

  GenotypeBuilder& Attribute(const std::string& key,
                             const std::string& value) {
    extended_attributes_[key] = value;
    return *this;
  }


 private:
  std::string sample_name_;
  std::vector<Allele*> alleles_;
  bool is_phased = false;
  int GQ_ = -1;
  int DP_ = -1;
  std::vector<int> AD_;
  std::vector<int> PL_;
  std::unordered_map<std::string, std::string> extended_attributes_;
};

} // easehts
} // ncic

#endif
