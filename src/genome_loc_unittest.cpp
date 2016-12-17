#include "easehts/genome_loc.h"
#include "easehts/unittest.h"
#include "easehts/reference_sequence.h"

#include <gtest/gtest.h>


using namespace ncic::easehts;

TEST(test_genome, normal) {
  const GenomeLoc& loc1 = GenomeLoc::kUnmapped;
  GenomeLoc loc2("loc2");
  GenomeLoc loc3("loc2", 3, 4, 2);
  GenomeLoc loc4("loc2", 4, 2, 1);
  GenomeLoc loc5("loc2", 4, 3, 1);
  GenomeLoc loc6("loc2", 4, 3, 2);
  EXPECT_STREQ(loc1.ToString().c_str(), "unmapped");
  EXPECT_TRUE(loc1 > loc3);
  EXPECT_TRUE(loc3 < loc4);
  EXPECT_TRUE(loc4 < loc5);
  EXPECT_TRUE(loc5 < loc6);

}

// TODO create reference test data
TEST(IntervalFileToList, IntervalUtils) {
  TEST_FILE("localtestdata/references/references/human_g1k_v37.fasta", b37);
  TEST_FILE("localtestdata/middle/mutect_panle_372.all.interval_list",
            interval_file);
  IndexedFastaSequenceFile reference(b37);
  GenomeLocParser parser(reference.GetSequenceDictionary());
  std::vector<GenomeLoc> genome_locs =
    IntervalUtils::IntervalFileToList(parser, interval_file, 0);
  EXPECT_EQ(genome_locs.size(), 16173);
  EXPECT_EQ(genome_locs[0], GenomeLoc("1", 0, 2488076, 2488196));
}
