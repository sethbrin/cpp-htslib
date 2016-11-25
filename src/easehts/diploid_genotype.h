//
// Created by zp on 11/25/16.
//

#ifndef EASEHTSLIB_DIPLOID_GENOTYPE_H_
#define EASEHTSLIB_DIPLOID_GENOTYPE_H_

#include "noncopyable.h"
#include "utils.h"
#include "base_utils.h"

#include <vector>

namespace ncic {
namespace easehts {

class DiploidGenotype {
 public:
  DiploidGenotype()
    :DiploidGenotype('N', 'N', -1) {}

  DiploidGenotype(char base1, char base2, int index)
    : base1_(base1),
    base2_(base2),
    index_(index) {}

  bool IsHomeRef(char r) const {
    return IsHom() && r == base1_;
  }

  bool IsHomeVar(char r) const {
    return IsHom() && r != base1_;
  }

  bool IsHetRef(char r) const {
    if (base1_ == r) return r != base2_;
    else return base2_ == r;
  }

  bool IsHom() const {
    return !IsHet();
  }

  bool IsHet() const {
    return base1_ != base2_;
  }

  int GetIndex() const {
    return index_;
  }

  /**
   * create a diploid genotype, given a character to make into a hom genotype
   *@param hom the character to turn into a hom genotype,
   * i.e. if it is A, then returned will be AA
   *@return the diploid genotype
   */
  static const DiploidGenotype& CreateHomGenotype(char hom) {
    int index = easehts::BaseUtils::SimpleBaseToBaseIndex(hom);
    ERROR_COND(index == -1,
               utils::StringFormatCStr("%c is not a valid base character", hom));
    return kGenotypes[kConversionMatrix[index][index]];
  }

  /**
   * create a diploid genotype, given 2 characters which may not necessarily be
   * ordered correctly
   * @param base1
   * @param base2
   * @return the diploid genotype
   */
  static const DiploidGenotype& CreateDiploidGenotype(char base1, char base2) {
    int index1 = easehts::BaseUtils::SimpleBaseToBaseIndex(base1);
    ERROR_COND(index1 == -1,
               utils::StringFormatCStr("%c is not a valid base character", base1));
    int index2 = easehts::BaseUtils::SimpleBaseToBaseIndex(base2);
    ERROR_COND(index2 == -1,
               utils::StringFormatCStr("%c is not a valid base character", base2));
    return kGenotypes[kConversionMatrix[index1][index2]];
  }

  /**
   * create a diploid genotype, given 2 characters which may not necessarily be
   * ordered correctly
   * @param base1
   * @param base2
   * @return the diploid genotype
   */
  static const DiploidGenotype& CreateDiploidGenotype(int base_index1, int base_index2) {
    ERROR_COND(base_index1 == -1,
               utils::StringFormatCStr("%d does no represent a valid character", base_index1));
    ERROR_COND(base_index2 == -1,
               utils::StringFormatCStr("%d does no represent a valid character", base_index2));
    return kGenotypes[kConversionMatrix[base_index1][base_index2]];
  }

  int Size() const {
    return kGenotypes.size();
  }

  const static std::vector<DiploidGenotype> kGenotypes;

 private:
  char base1_;
  char base2_;

  // the index in the kConversionMatrix, start from 0
  int index_;

  const static int kConversionMatrix[4][4];
};


} // easehts
} // ncic

#endif
