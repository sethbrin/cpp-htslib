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

void unordered_map_normal(
    const std::vector<gatk::PileupElement*>& elements) {
  int size = elements.size();
  // store the read_name=>read_idx
  std::unordered_map<std::string, int> filter_map;

  std::vector<bool> elements_to_keep(size, false);

  for (int index=0; index<size; index++) {
    const std::string& read_name = elements[index]->GetRead()->GetQueryName();

    if (filter_map.find(read_name) == filter_map.end()) {
      filter_map[read_name] = index;
      elements_to_keep[index] = true;
    } else {
      int existing_index = filter_map[read_name];
      const gatk::PileupElement& cur_pileup = *elements[index];
      const gatk::PileupElement& existing_pileup = *elements[existing_index];

      // if the read disagree at this position
      if (existing_pileup.GetBase() != cur_pileup.GetBase()) {
        // keep the mismatching one
        if (cur_pileup.GetBase() != 'A') {
          elements_to_keep[existing_index] = false;
          filter_map[read_name] = index;
          elements_to_keep[index] = true;
        }
      } else {
        // otherwise, keep the element with the higher qualitu score
        if (existing_pileup.GetQual() < cur_pileup.GetQual()) {
          elements_to_keep[existing_index] = false;
          filter_map[read_name] = index;
          elements_to_keep[index] = true;
        }
      }
    }
  }

}

void map_normal(
    const std::vector<gatk::PileupElement*>& elements) {
  int size = elements.size();
  // store the read_name=>read_idx
  std::map<std::string, int> filter_map;

  std::vector<bool> elements_to_keep(size, false);

  for (int index=0; index<size; index++) {
    const std::string& read_name = elements[index]->GetRead()->GetQueryName();

    if (filter_map.find(read_name) == filter_map.end()) {
      filter_map[read_name] = index;
      elements_to_keep[index] = true;
    } else {
      int existing_index = filter_map[read_name];
      const gatk::PileupElement& cur_pileup = *elements[index];
      const gatk::PileupElement& existing_pileup = *elements[existing_index];

      // if the read disagree at this position
      if (existing_pileup.GetBase() != cur_pileup.GetBase()) {
        // keep the mismatching one
        if (cur_pileup.GetBase() != 'A') {
          elements_to_keep[existing_index] = false;
          filter_map[read_name] = index;
          elements_to_keep[index] = true;
        }
      } else {
        // otherwise, keep the element with the higher qualitu score
        if (existing_pileup.GetQual() < cur_pileup.GetQual()) {
          elements_to_keep[existing_index] = false;
          filter_map[read_name] = index;
          elements_to_keep[index] = true;
        }
      }
    }
  }

}

struct mystring_hash {
  size_t operator()(const std::string& s) const {
    unsigned long __h = 0;
    for (unsigned i = 0;i < s.size();++i)
      __h = 31 * __h + s[i];
    return size_t(__h);
  }
};

void unordered_map_user_hash(
    const std::vector<gatk::PileupElement*>& elements) {
  int size = elements.size();
  // store the read_name=>read_idx
  std::unordered_map<std::string, int, mystring_hash> filter_map;

  std::vector<bool> elements_to_keep(size, false);

  for (int index=0; index<size; index++) {
    const std::string& read_name = elements[index]->GetRead()->GetQueryName();

    if (filter_map.find(read_name) == filter_map.end()) {
      filter_map[read_name] = index;
      elements_to_keep[index] = true;
    } else {
      int existing_index = filter_map[read_name];
      const gatk::PileupElement& cur_pileup = *elements[index];
      const gatk::PileupElement& existing_pileup = *elements[existing_index];

      // if the read disagree at this position
      if (existing_pileup.GetBase() != cur_pileup.GetBase()) {
        // keep the mismatching one
        if (cur_pileup.GetBase() != 'A') {
          elements_to_keep[existing_index] = false;
          filter_map[read_name] = index;
          elements_to_keep[index] = true;
        }
      } else {
        // otherwise, keep the element with the higher qualitu score
        if (existing_pileup.GetQual() < cur_pileup.GetQual()) {
          elements_to_keep[existing_index] = false;
          filter_map[read_name] = index;
          elements_to_keep[index] = true;
        }
      }
    }
  }

}

struct record_hash {
  size_t operator()(SAMBAMRecord* read) const {
    return read->HashCode();
  }
};

struct record_equal {
  bool operator()(SAMBAMRecord* lhs,
                  SAMBAMRecord* rhs) const {
    return lhs->GetQueryName() == rhs->GetQueryName();
  }
};

void unordered_map_hash_read(
    const std::vector<gatk::PileupElement*>& elements) {
  int size = elements.size();
  // store the read=>read_idx
  std::unordered_map<SAMBAMRecord*, int,
    record_hash, record_equal> filter_map;

  std::vector<bool> elements_to_keep(size, false);

  for (int index=0; index<size; index++) {
    SAMBAMRecord* read = elements[index]->GetRead();

    if (filter_map.find(read) == filter_map.end()) {
      filter_map[read] = index;
      elements_to_keep[index] = true;
    } else {
      int existing_index = filter_map[read];
      const gatk::PileupElement& cur_pileup = *elements[index];
      const gatk::PileupElement& existing_pileup = *elements[existing_index];

      // if the read disagree at this position
      if (existing_pileup.GetBase() != cur_pileup.GetBase()) {
        // keep the mismatching one
        if (cur_pileup.GetBase() != 'A') {
          elements_to_keep[existing_index] = false;
          filter_map[read] = index;
          elements_to_keep[index] = true;
        }
      } else {
        // otherwise, keep the element with the higher qualitu score
        if (existing_pileup.GetQual() < cur_pileup.GetQual()) {
          elements_to_keep[existing_index] = false;
          filter_map[read] = index;
          elements_to_keep[index] = true;
        }
      }
    }
  }

}


GenomeLoc interval("X", 22, 39911326, 39911335);

static void BM_unordered_map_normal(benchmark::State& state) {
  TEST_FILE("MG225_normal_sorted_X.bam", filename);
  BAMIndexReader reader(filename);
  std::vector<gatk::PileupElement*> elements;
  int pos = state.range(0);
  gatk::GATKPileupTraverse traverse(&reader, interval);

  while (state.KeepRunning()) {
    state.PauseTiming();
    PreparePileup(traverse, elements, pos, interval.GetStop());
    state.ResumeTiming();
    unordered_map_normal(elements);
  }
}

static void BM_map_normal(benchmark::State& state) {
  TEST_FILE("MG225_normal_sorted_X.bam", filename);
  BAMIndexReader reader(filename);
  std::vector<gatk::PileupElement*> elements;
  int pos = state.range(0);
  gatk::GATKPileupTraverse traverse(&reader, interval);

  while (state.KeepRunning()) {
    state.PauseTiming();
    PreparePileup(traverse, elements, pos, interval.GetStop());
    state.ResumeTiming();
    map_normal(elements);
  }
}

static void BM_unordered_map_user_hash(benchmark::State& state) {
  TEST_FILE("MG225_normal_sorted_X.bam", filename);
  BAMIndexReader reader(filename);
  std::vector<gatk::PileupElement*> elements;
  int pos = state.range(0);
  gatk::GATKPileupTraverse traverse(&reader, interval);

  while (state.KeepRunning()) {
    state.PauseTiming();
    PreparePileup(traverse, elements, pos, interval.GetStop());
    state.ResumeTiming();
    unordered_map_user_hash(elements);
  }
}

static void BM_unordered_map_hash_read(benchmark::State& state) {
  TEST_FILE("MG225_normal_sorted_X.bam", filename);
  BAMIndexReader reader(filename);
  std::vector<gatk::PileupElement*> elements;
  int pos = state.range(0);
  gatk::GATKPileupTraverse traverse(&reader, interval);

  while (state.KeepRunning()) {
    state.PauseTiming();
    PreparePileup(traverse, elements, pos, interval.GetStop());
    state.ResumeTiming();
    unordered_map_hash_read(elements);
  }
}



// Register the function as a benchmark
BENCHMARK(BM_unordered_map_normal)->DenseRange(
    interval.GetStart(),
    interval.GetStop(), 1);

BENCHMARK(BM_map_normal)->DenseRange(
    interval.GetStart(),
    interval.GetStop(), 1);

BENCHMARK(BM_unordered_map_user_hash)->DenseRange(
    interval.GetStart(),
    interval.GetStop(), 1);

BENCHMARK(BM_unordered_map_hash_read)->DenseRange(
    interval.GetStart(),
    interval.GetStop(), 1);

BENCHMARK_MAIN();
