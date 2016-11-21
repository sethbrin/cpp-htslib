//
// Created by zp on 11/19/16.
//

#include "easehts/easehts.h"
#include "easehts/unittest.h"

#include <gtest/gtest.h>
#include <stdlib.h>
#include <memory>

using namespace ncic::easehts;

samFile* GetSamFile() {
  TEST_FILE("uncompressed.sam", filename);
  samFile* in = sam_open(filename, "r");
  if (in == NULL) {
    fprintf(stderr, "Error opening \"%s\"\n", filename);
    exit(1);
  }
  return in;
}

TEST(Constructor1, SAMBAMTextReader) {
  SAMBAMTextReader reader(GetSamFile());
  SAMBAMRecord record;
  int count = 0;
  while (reader.HasNext(&record)) {
    count++;
  }

  EXPECT_EQ(count, 10);

}

TEST(MoveConstructor, SAMBAMTextReader) {
  std::vector<SAMBAMTextReader> readers;
  SAMBAMTextReader reader = SAMBAMTextReader(GetSamFile());
  readers.push_back(std::move(reader));
  SAMBAMRecord record;
  int count = 0;
  while (readers[0].HasNext(&record)) {
    count++;
  }

  EXPECT_EQ(count, 10);

}

TEST(Constructor2, SAMBAMTextReader) {
  TEST_FILE("uncompressed.sam", filename);
  SAMBAMTextReader reader(filename);
  SAMBAMRecord record;
  int count = 0;
  while (reader.HasNext(&record)) {
    count++;
  }

  EXPECT_EQ(count, 10);
}

TEST(AddRegion, BAMIndexReader) {
  TEST_FILE("index_test.bam", filename);
  BAMIndexReader reader(filename);
  reader.SetRegion("chr1");
  SAMBAMRecord record;
  int count = 0;
  while (reader.HasNext(&record)) {
    count++;
  }

  EXPECT_EQ(count, 885);
}

TEST(AddRegion2, BAMIndexReader) {
  TEST_FILE("index_test.bam", filename);
  BAMIndexReader reader(filename);
  reader.SetRegion("chr1");
  SAMBAMRecord record;
  int count = 0;
  while (reader.HasNext(&record)) {
    if (count == 10) break;
    count++;
  }
  EXPECT_EQ(count, 10);

  reader.SetRegion("chr2");
  count = 0;
  while (reader.HasNext(&record)) {
    count++;
  }
  EXPECT_EQ(count, 837);
}

TEST(WithoutSetRegion, BAMIndexReader) {
  TEST_FILE("index_test.bam", filename);
  BAMIndexReader reader(filename);
  int count = 0;
  SAMBAMRecord record;
  while (reader.HasNext(&record)) {
    count++;
  }
  EXPECT_EQ(count, 10000);
}

TEST(AddRegion, BAMIndexBatchReader) {
  TEST_FILE("index_test.bam", filename);
  BAMIndexBatchReader reader(filename);
  reader.AddRegion("chr2");
  reader.AddRegion("chr1");
  SAMBAMRecord record;
  int count = 0;
  while (reader.HasNext(&record)) {
    count++;
  }

  EXPECT_EQ(count, 1722);
}

TEST(BAMIndexReaderTest, BamRecord) {
  TEST_FILE("index_test.bam", filename);
  BAMIndexBatchReader reader(filename);
  reader.AddRegion("chr2:236164650-236164650");
  SAMBAMRecord record;
  if (reader.HasNext(&record)) {
    EXPECT_STREQ(record.GetQueryName(), "5934311");
    EXPECT_STREQ(record.GetSequence().c_str(), "CTTGAACTCCTGGACTCAGGTGATGCATCCGCCTCAGCCTCCCAAACTGTT");
  }
}
