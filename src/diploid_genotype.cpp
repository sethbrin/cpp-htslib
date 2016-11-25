#include "easehts/diploid_genotype.h"


namespace ncic {
namespace easehts {

const std::vector<DiploidGenotype> DiploidGenotype::kGenotypes = {
  DiploidGenotype('A', 'A', 0),
  DiploidGenotype('A', 'C', 1),
  DiploidGenotype('C', 'C', 2),
  DiploidGenotype('A', 'G', 3),
  DiploidGenotype('C', 'G', 4),
  DiploidGenotype('G', 'G', 5),
  DiploidGenotype('A', 'T', 6),
  DiploidGenotype('C', 'T', 7),
  DiploidGenotype('G', 'T', 8),
  DiploidGenotype('T', 'T', 9)
};

const int DiploidGenotype::kConversionMatrix[4][4] = {
  {0, 1, 3, 6},
  {1, 2, 4, 7},
  {3, 4, 5, 8},
  {6, 7, 8, 9},
};


} // easehts
} // ncic
