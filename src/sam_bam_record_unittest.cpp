//
// Created by zp on 9/5/16.
//
#include "easehts/easehts.h"
#include "easehts/unittest.h"

#include <gtest/gtest.h>
#include <memory>

using namespace ncic::easehts;

TEST(SAMBAMRecord, recordCopy) {
  TEST_FILE("uncompressed.sam", filename);
  SAMBAMTextReader reader(filename);

  SAMBAMRecord record;
  SAMBAMRecord record2(false);

  if (reader.HasNext(&record)) {
    EXPECT_EQ(record2.GetRawRecord(), reinterpret_cast<void*>(NULL));
    record2 = record.Copy();
    EXPECT_STREQ(record2.GetSequence().c_str(), "CAACAGAAGC");
    EXPECT_STREQ(record.GetSequence().c_str(), "CAACAGAAGC");

    // copy
    SAMBAMRecord record3 = record.Copy();
    EXPECT_STREQ(record3.GetSequence().c_str(), "CAACAGAAGC");

    // move
    SAMBAMRecord record4(std::move(record));
    EXPECT_STREQ(record3.GetSequence().c_str(), "CAACAGAAGC");
  }
}

TEST(SAMFlag, value) {
  EXPECT_EQ(SAMFlag::READ_PAIRED, 0x1);
  EXPECT_EQ(SAMFlag::NOT_PRIMARY_ALIGNMENT, 0x100);
}
