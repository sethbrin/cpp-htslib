#include "easehts/genotype.h"

namespace ncic {
namespace easehts {

const std::string kNoCallString = ".";
const std::string kSpanDelString = "*";

const std::shared_ptr<Allele> Allele::kRefA = std::make_shared<Allele>("A", true);
const std::shared_ptr<Allele> Allele::kAltA = std::make_shared<Allele>("A", false);
const std::shared_ptr<Allele> Allele::kRefC = std::make_shared<Allele>("C", true);
const std::shared_ptr<Allele> Allele::kAltC = std::make_shared<Allele>("C", false);
const std::shared_ptr<Allele> Allele::kRefG = std::make_shared<Allele>("G", true);
const std::shared_ptr<Allele> Allele::kAltG = std::make_shared<Allele>("G", false);
const std::shared_ptr<Allele> Allele::kRefT = std::make_shared<Allele>("T", true);
const std::shared_ptr<Allele> Allele::kAltT = std::make_shared<Allele>("T", false);
const std::shared_ptr<Allele> Allele::kRefN = std::make_shared<Allele>("N", true);
const std::shared_ptr<Allele> Allele::kAltN = std::make_shared<Allele>("N", false);
const std::shared_ptr<Allele> Allele::kNoCall = std::make_shared<Allele>(".", false);
const std::shared_ptr<Allele> Allele::kSpanDel = std::make_shared<Allele>("*", false);

} // easehts
} // ncic
