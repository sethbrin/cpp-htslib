//
// Created by zp on 8/25/16.
//

#include "easehts.h"

#include <gtest/gtest.h>
#include <stdlib.h>
#include <memory>

samFile* GetSamFile() {
  char filename[] = "/home/zp/work/rnaseq/cpp-htslib/test/uncompressed.sam";
  samFile* in = sam_open(filename, "r");
  if (in == NULL) {
    fprintf(stderr, "Error opening \"%s\"\n", filename);
    exit(1);
  }
  return in;
}

ncic::easehts::SAMBAMNormalReader reader(GetSamFile());

TEST(SAMBAMNormalReaderTest, Sam) {
  ncic::easehts::SAMBAMRecord record;
  int count = 0;
  while (reader.HasNext(&record)) {
    count++;
  }

  EXPECT_EQ(count, 10);
}

samFile* GetBamFile() {
  char filename[] = "/home/zp/work/rnaseq/cpp-htslib/test/compressed.bam";
  samFile* in = sam_open(filename, "r");
  if (in == NULL) {
    fprintf(stderr, "Error opening \"%s\"\n", filename);
    exit(1);
  }
  return in;
}

ncic::easehts::SAMBAMNormalReader reader2(GetBamFile());
TEST(SAMBAMNormalReaderTest, bam) {
  ncic::easehts::SAMBAMRecord record;
  int count = 0;
  while (reader2.HasNext(&record)) {
    count++;
  }

  EXPECT_EQ(count, 10);
}

samFile* GetBamIndexFile() {
  char filename[] = "/home/zp/work/rnaseq/cpp-htslib/test/index_test.bam";
  samFile* in = sam_open(filename, "r");
  if (in == NULL) {
    fprintf(stderr, "Error opening \"%s\"\n", filename);
    exit(1);
  }
  return in;
}

std::shared_ptr<ncic::easehts::BAMIndex> bam_index =
    std::make_shared<ncic::easehts::BAMIndex>("/home/zp/work/rnaseq/cpp-htslib/test/index_test.bam");

ncic::easehts::BAMIndexReader reader_index(GetBamIndexFile(), bam_index);
TEST(BAMIndexReaderTest, BamIndex) {
  reader_index.AddRegion("chr1");
  ncic::easehts::SAMBAMRecord record;
  int count = 0;
  while (reader_index.HasNext(&record)) {
    count++;
  }

  EXPECT_EQ(count, 885);
}

TEST(BAMIndexReaderTest, BamIndex2) {
  reader_index.AddRegion("chr2");
  reader_index.AddRegion("chr1");
  ncic::easehts::SAMBAMRecord record;
  int count = 0;
  while (reader_index.HasNext(&record)) {
    count++;
  }

  EXPECT_EQ(count, 1722);
}

TEST(BAMIndexReaderTest, BamIndex3) {
  reader_index.AddRegion("chr2:236164650-242159443");
  ncic::easehts::SAMBAMRecord record;
  int count = 0;
  while (reader_index.HasNext(&record)) {
    count++;
  }

  EXPECT_EQ(count, 27);
}

TEST(BAMIndexReaderTest, BamIndex4) {
  reader_index.AddRegion("chr2:236164650-242159443");
  ncic::easehts::SAMBAMRecord record;
  int count = 0;
  while (reader_index.HasNext(&record)) {
    count++;
  }

  EXPECT_EQ(count, 27);
}

TEST(BAMIndexReaderTest, BamRecord) {
  reader_index.AddRegion("chr2:236164650-236164650");
  ncic::easehts::SAMBAMRecord record;
  if (reader_index.HasNext(&record)) {
    EXPECT_STREQ(record.GetQueryName(), "5934311");
    EXPECT_STREQ(record.GetSequence().c_str(), "CTTGAACTCCTGGACTCAGGTGATGCATCCGCCTCAGCCTCCCAAACTGTT");
  }
}
