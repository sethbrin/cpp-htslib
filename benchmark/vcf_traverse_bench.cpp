#include <benchmark/benchmark.h>
#include <vector>
#include <unordered_map>
#include <map>

#include <easehts/vcf.h>
#include <easehts/unittest.h>
#include <easehts/genome_loc.h>

using namespace ncic::easehts;

static void vcf_bench() {
  TEST_FILE(
      "localtestdata/references/vcf/dbsnp_138.b37.compressed.vcf",
      filename);
  VCFIndexReader reader(filename);
  reader.SetRegion(22, 39931576, 39934456);

  VCFTraverse traverse(&reader);

  for (int pos = 39931576; pos <= 39934456; pos++) {
    traverse.SeekFroward(GenomeLoc("X", 22, pos, pos));
  }
}

static void BM_vcf(benchmark::State& state) {
  while (state.KeepRunning()) {
    vcf_bench();
  }
}

BENCHMARK(BM_vcf);

BENCHMARK_MAIN();

