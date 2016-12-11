#include "easehts/unittest.h"
#include "easehts/pileup.h"
#include "easehts/sam_bam_reader.h"
#include "easehts/sam_bam_record.h"

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
    ReadBackedRawPileup plp = traverse.Next();
    if (count == 0) {
      EXPECT_EQ(plp.Size(), 5);
    }
    count ++;
    //printf("contig:%d pos:%d count:%d\n", plp.GetContigId(),
    //       plp.GetPos(), plp.Size());
  }

  EXPECT_EQ(traverse.Next().Size(), 0);
  EXPECT_EQ(count, 4101);
}

TEST(MutectBam, NormalPileupTraverse) {
  TEST_FILE("MG225_normal_sorted_X.bam", filename);
  BAMIndexReader reader(filename);
  reader.SetRegion(22, 15482483-3000, 15482603+3000);
  TraverseData data;
  data.reader = &reader;

  PileupTraverse traverse(fun, &data);
  while (traverse.HasNext()) {
    ReadBackedRawPileup plp = traverse.Next();
    if (plp.GetPos() == 15482483) {
      EXPECT_EQ(plp.Size(), 158);
      break;
    }
  }
}

TEST(MutectBam, TumorPileupTraverse) {
  TEST_FILE("MG225_tumor_sorted_X.bam", filename);
  BAMIndexReader reader(filename);
  reader.SetRegion(22, 15482483-3000, 15482603+3000);
  TraverseData data;
  data.reader = &reader;

  PileupTraverse traverse(fun, &data);
  while (traverse.HasNext()) {
    ReadBackedRawPileup plp = traverse.Next();
    if (plp.GetPos() == 15482483) {
      EXPECT_EQ(plp.Size(), 266);
      break;
    }
  }

}

TEST(GATKPileupElement, Destruction) {
  TEST_FILE("MG225_tumor_sorted_X.bam", filename);
  BAMIndexReader reader(filename);
  reader.SetRegion(22, 15482483-3000, 15482603+3000);

  SAMBAMRecord* record = new SAMBAMRecord();
  reader.HasNext(record);

  GATKPileupElement* element = new GATKPileupElement(record, 0, 0, false, 1);
  EXPECT_NE(element->GetRead()->GetRawRecord(), nullptr);
  delete element;
}
