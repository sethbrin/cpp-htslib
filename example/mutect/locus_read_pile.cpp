#include "locus_read_pile.h"

namespace ncic {
namespace mutect {


void LocusReadPile::AddPileupElement(const easehts::ReadBackedPileup& read_backed_pileup) {
  for (size_t i=0; i < read_backed_pileup.Size(); i++) {
    elements_.emplace_back(read_backed_pileup[i]);
  }
}

} // mutect
} // ncic
