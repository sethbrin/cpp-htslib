#include "unittest.h"
#include "pileup.h"
#include "sam_bam_reader.h"
#include "sam_bam_record.h"

#include <gtest/gtest.h>

using namespace ncic::easehts;

typedef struct TraverseData {
  BAMIndexReader* reader;
} TraverseData;

int fun(void *data, bam1_t *b) {
  TraverseData& traverse_data = *(TraverseData*)data;
  while (traverse_data.reader->HasNext(b)) {
    return 0;
  }
  return -1;
}

TEST(Next, PileupTraverse) {
  TEST_FILE("pileup/mpileup.1.bam", filename);
  BAMIndexReader reader(filename);
  TraverseData data;
  data.reader = &reader;

  PileupTraverse traverse(fun, &data);
  int count = 0;
  while (traverse.HasNext()) {
    ReadBackedPileup plp = traverse.Next();
    if (count == 0) {
      EXPECT_EQ(plp.Size(), 5);
    }
    count ++;
    //printf("contig:%d pos:%d count:%d\n", plp.GetContigId(),
    //       plp.GetPos(), plp.Size());
  }
  EXPECT_EQ(count, 4101);
}
