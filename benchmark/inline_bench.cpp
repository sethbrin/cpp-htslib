#include <benchmark/benchmark.h>
#include <vector>
#include <unordered_map>
#include <map>

#include <easehts/gatk/pileup.h>
#include <easehts/sam_bam_reader.h>
#include <easehts/sam_bam_record.h>
#include <easehts/unittest.h>
#include <easehts/genome_loc.h>

using namespace ncic::easehts;

using PileupFilterFun = std::function<bool (gatk::PileupElement*)>;
/**
 * The filters to filter gatk::PileupElement
 *
 * when the PileupFilterFun returns true, the element remains,
 * otherwise removed
 */
class PileupFilter : public NonCopyable {
 public:
  PileupFilter() = delete;
  // ignore the element which is deletion
  static bool IsNotDeletion(gatk::PileupElement* element) {
    return !element->IsDeletion();
  }

  /**
   * when the quality >= min_base_quality and
   * map quality >= min_map_quality returns true
   */
  static bool IsBaseAndMappingQualityLarge(
      gatk::PileupElement* element,
      int min_base_quality, int min_map_quality) {
    bool a =  element->GetRead()->GetMapQuality() >= min_map_quality;
    bool b = (element->IsDeletion() || element->GetQual() >= min_base_quality);
    return element->GetRead()->GetMapQuality() >= min_map_quality &&
      (element->IsDeletion() || element->GetQual() >= min_base_quality);
  }

  static bool IsBaseQualityLarge(
      gatk::PileupElement* element, int min_base_quality) {
    return (element->IsDeletion() || element->GetQual() >= min_base_quality);
  }

  static bool IsMappingQualityLarger(
      gatk::PileupElement* element,
      int min_map_quality) {
    return element->GetRead()->GetMapQuality() >= min_map_quality;
  }

  static bool IsMappingQualityLargerThanZero(gatk::PileupElement* element) {
    return element->GetRead()->GetMapQuality() > 0;
  }


  static bool IsPositiveStrand(gatk::PileupElement* element) {
    return !element->GetRead()->GetReadNegativeStrandFlag();
  }

  static bool IsNegativeStrand(gatk::PileupElement* element) {
    return element->GetRead()->GetReadNegativeStrandFlag();
  }
};

static inline void GetPileupByFilter(
    const std::vector<gatk::PileupElement*>& elements,
    gatk::ReadBackedPileup* pPileup,
    PileupFilterFun pred) {
  assert(pPileup->Size() == 0);
  int size = elements.size();
  for (int i = 0; i < size; i++) {
    if (pred(elements[i])) {
      pPileup->AddElement(elements[i]);
    }
  }
}

#ifdef GET_PILEUP_BY_FILTER
#undef GET_PILEUP_BY_FILTER
#endif

#define GET_PILEUP_BY_FILTER(elements, pileup, pred) \
  assert(pileup.Size() == 0); \
  int size = elements.size(); \
  for (int i = 0; i < elements.size(); i++) { \
    if (pred) { \
      pileup.AddElement(elements[i]); \
    } \
  }

bool PreparePileup(
    gatk::GATKPileupTraverse& traverse,
    std::vector<gatk::PileupElement*>& elements,
    int pos,
    int end_pos) {
  elements.clear();
  while (traverse.HasNext()) {
    if (traverse.CurrentPileup().GetPos() == pos) {
      elements = traverse.CurrentPileup().GetElements();
      return true;
    }
    if (traverse.CurrentPileup().GetPos() > end_pos) return false;
  }
  return false;
}

static void normal(const std::vector<gatk::PileupElement*>& elements) {
  gatk::ReadBackedPileup newpile;

  GetPileupByFilter(elements, &newpile, std::bind(
          PileupFilter::IsBaseAndMappingQualityLarge,
          std::placeholders::_1, -1, -1));
}

static void use_inline(const std::vector<gatk::PileupElement*>& elements) {
  gatk::ReadBackedPileup newpile;

  int size = elements.size();
  for (int i = 0; i < size; i++) {
    if (PileupFilter::IsBaseAndMappingQualityLarge(elements[i], -1, -1)) {
      newpile.AddElement(elements[i]);
    }
  }
}

static void use_define(const std::vector<gatk::PileupElement*>& elements) {
  gatk::ReadBackedPileup newpile;

  GET_PILEUP_BY_FILTER(
      elements, newpile,
      PileupFilter::IsBaseAndMappingQualityLarge(elements[i], -1, -1));

}


GenomeLoc interval("X", 22, 39911326, 39911335);

static void BM_normal(benchmark::State& state) {
  TEST_FILE("MG225_normal_sorted_X.bam", filename);
  BAMIndexReader reader(filename);
  std::vector<gatk::PileupElement*> elements;
  int pos = state.range(0);
  gatk::GATKPileupTraverse traverse(&reader, interval);

  while (state.KeepRunning()) {
    state.PauseTiming();
    PreparePileup(traverse, elements, pos, interval.GetStop());
    state.ResumeTiming();
    normal(elements);
  }
}

static void BM_inline(benchmark::State& state) {
  TEST_FILE("MG225_normal_sorted_X.bam", filename);
  BAMIndexReader reader(filename);
  std::vector<gatk::PileupElement*> elements;
  int pos = state.range(0);
  gatk::GATKPileupTraverse traverse(&reader, interval);

  while (state.KeepRunning()) {
    state.PauseTiming();
    PreparePileup(traverse, elements, pos, interval.GetStop());
    state.ResumeTiming();
    use_inline(elements);
  }
}

static void BM_define(benchmark::State& state) {
  TEST_FILE("MG225_normal_sorted_X.bam", filename);
  BAMIndexReader reader(filename);
  std::vector<gatk::PileupElement*> elements;
  int pos = state.range(0);
  gatk::GATKPileupTraverse traverse(&reader, interval);

  while (state.KeepRunning()) {
    state.PauseTiming();
    PreparePileup(traverse, elements, pos, interval.GetStop());
    state.ResumeTiming();
    use_define(elements);
  }
}


// Register the function as a benchmark
BENCHMARK(BM_normal)->DenseRange(
    interval.GetStart(),
    interval.GetStop(), 1);

BENCHMARK(BM_inline)->DenseRange(
    interval.GetStart(),
    interval.GetStop(), 1);

BENCHMARK(BM_define)->DenseRange(
    interval.GetStart(),
    interval.GetStop(), 1);

BENCHMARK_MAIN();
