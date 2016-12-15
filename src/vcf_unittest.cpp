#include "easehts/vcf.h"
#include "easehts/unittest.h"

#include <gtest/gtest.h>

using namespace ncic::easehts;

TEST(VCFReader, VCFHeader) {
  TEST_FILE("vcf/exampleDBSNP.vcf", filename);
  VCFReader reader(filename);

  const VCFHeader& header = reader.GetHeader();

  EXPECT_STREQ(header.GetVersion().c_str(), "VCFv4.1");
  EXPECT_EQ(header.GetSamplesCnt(), 0);

  int cnt = 0;
  while (reader.HasNext()) {
    VariantContext* record = reader.Next();
    cnt ++;
    delete record;
  }
  EXPECT_EQ(cnt, 217);
}

