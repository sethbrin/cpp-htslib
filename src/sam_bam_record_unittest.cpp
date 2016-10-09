//
// Created by zp on 9/5/16.
//
#include "easehts.h"
#include "unittest.h"

#include <gtest/gtest.h>
#include <memory>

samFile* GetSamFile() {
  TEST_FILE("uncompressed.sam", filename);
  samFile* in = sam_open(filename, "r");
  if (in == NULL) {
    fprintf(stderr, "Error opening \"%s\"\n", filename);
    exit(1);
  }
  return in;
}

ncic::easehts::SAMBAMNormalReader reader(GetSamFile());

TEST(SAMBAMRecord, recordCopy) {
  ncic::easehts::SAMBAMRecord record;
  ncic::easehts::SAMBAMRecord record2(false);

  if (reader.HasNext(&record)) {
    EXPECT_EQ(record2.GetRawRecord(), reinterpret_cast<void*>(NULL));
    record2 = record;
    EXPECT_STREQ(record2.GetSequence().c_str(), "CAACAGAAGC");
    EXPECT_STREQ(record.GetSequence().c_str(), "CAACAGAAGC");

    // copy
    ncic::easehts::SAMBAMRecord record3(record);
    EXPECT_STREQ(record3.GetSequence().c_str(), "CAACAGAAGC");

    // move
    ncic::easehts::SAMBAMRecord record4(std::move(record));
    EXPECT_STREQ(record3.GetSequence().c_str(), "CAACAGAAGC");
  }
}
