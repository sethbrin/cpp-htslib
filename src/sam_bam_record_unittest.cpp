//
// Created by zp on 9/5/16.
//
#include "easehts/easehts.h"
#include "easehts/unittest.h"

#include <gtest/gtest.h>
#include <memory>
#include <unordered_map>

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

TEST(SAMBAMRecord, CigarElement) {
  CigarElement element('M', 100);
  EXPECT_EQ(element.ToString(), "100M");
}

TEST(SAMBAMRecord, GetCigarString) {
  TEST_FILE("uncompressed.sam", filename);
  SAMBAMTextReader reader(filename);

  SAMBAMRecord record;

  if (reader.HasNext(&record)) {
    EXPECT_EQ(record.GetCigarString(), "8M2I");
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

TEST(SAMBAMRecord, unordered_map) {
  TEST_FILE("uncompressed.sam", filename);
  SAMBAMTextReader reader(filename);

  SAMBAMRecord record;

  reader.HasNext(&record);
  SAMBAMRecord record2 = record.Copy();
  std::unordered_map<SAMBAMRecord*, int,
    record_hash, record_equal> map;
  EXPECT_EQ(map.find(&record), map.end());
  map[&record]  = 1;
  EXPECT_NE(map.find(&record), map.end());
  EXPECT_EQ(map[&record], 1);
  EXPECT_EQ(map[&record2], 1);
  EXPECT_NE(map.find(&record2), map.end());
  SAMBAMRecord record3;
  reader.HasNext(&record3);
  EXPECT_NE(map.find(&record3), map.end());
  EXPECT_EQ(map[&record3], 1);
  map[&record3] = 3;
  EXPECT_NE(map.find(&record3), map.end());
}
