#include "easehts/vcf.h"
#include "easehts/unittest.h"

#include <thread>
#include <vector>

#include <gtest/gtest.h>

using namespace ncic::easehts;

TEST(VCFTextReader, normal) {
  TEST_FILE("vcf/exampleDBSNP.vcf", filename);
  VCFTextReader reader(filename);

  int cnt = 0;
  while (reader.HasNext()) {
    VariantContext* record = reader.Next();
    cnt ++;
    delete record;
  }
  EXPECT_EQ(cnt, 217);
}

TEST(VCFIndexReader, normal) {
  TEST_FILE("vcf/exampleDBSNP.vcf", filename);
  VCFIndexReader reader(filename);

  int cnt = 0;
  while (reader.HasNext()) {
    VariantContext* record = reader.Next();
    cnt ++;
    delete record;
  }
  EXPECT_EQ(cnt, 217);
}

TEST(VCFIndexReader, SetRegion) {
  TEST_FILE("vcf/exampleDBSNP.compress.vcf", filename);
  VCFIndexReader reader(filename);
  reader.SetRegion("chr1:10440-10519");

  int cnt = 0;
  while (reader.HasNext()) {
    VariantContext* record = reader.Next();
    cnt ++;
    delete record;
  }
  EXPECT_EQ(cnt, 5);

  reader.SetRegion(0, 10439, 10519);
  cnt = 0;
  while (reader.HasNext()) {
    VariantContext* record = reader.Next();
    cnt ++;
    delete record;
  }
  EXPECT_EQ(cnt, 5);

}



TEST(VCFHeader, normal) {
  TEST_FILE("vcf/exampleDBSNP.vcf", filename);
  VCFTextReader reader(filename);

  const VCFHeader& header = reader.GetHeader();

  EXPECT_STREQ(header.GetVersion().c_str(), "VCFv4.1");
  EXPECT_EQ(header.GetSamplesCnt(), 0);


  EXPECT_EQ(header.GetContigId("chr1"), 0);
  EXPECT_STREQ(header.GetContig(0).c_str(), "chr1");
  EXPECT_EQ(header.GetContigLength(0), 249250621);
  EXPECT_EQ(header.GetContigLength("chr1"), 249250621);
  EXPECT_STREQ(header.CreateOverEntireContig(0).ToString().c_str(),
            ":0-249250620");
}

TEST(VCFTraverse, normal) {
  TEST_FILE("vcf/exampleDBSNP.vcf", filename);
  VCFTextReader reader(filename);

  VCFTraverse traverse(&reader);

  traverse.SeekFroward(GenomeLoc("chr1", 0, 10438, 10438));
  std::vector<VariantContext*> records;
  traverse.GetRecords(&records);
  EXPECT_EQ(records.size(), 2);

  traverse.SeekFroward(GenomeLoc("chr1", 0, 10439, 10439));
  traverse.GetRecords(&records);
  EXPECT_EQ(records.size(), 1);

  traverse.SeekFroward(GenomeLoc("chr1", 0, 10440, 10440));
  traverse.GetRecords(&records);
  EXPECT_EQ(records.size(), 0);


  traverse.SeekFroward(GenomeLoc("chr1", 0, 10440, 10519));
  traverse.GetRecords(&records);
  EXPECT_EQ(records.size(), 2);
}

TEST(VCFTraverse, VCFIndexReader) {
  TEST_FILE("vcf/exampleDBSNP.compress.vcf", filename);
  VCFIndexReader reader(filename);
  reader.SetRegion(0, 10438, 10519);

  VCFTraverse traverse(&reader);

  traverse.SeekFroward(GenomeLoc("chr1", 0, 10438, 10438));
  std::vector<VariantContext*> records;
  traverse.GetRecords(&records);
  EXPECT_EQ(records.size(), 2);

  traverse.SeekFroward(GenomeLoc("chr1", 0, 10439, 10439));
  traverse.GetRecords(&records);
  EXPECT_EQ(records.size(), 1);

  traverse.SeekFroward(GenomeLoc("chr1", 0, 10440, 10440));
  traverse.GetRecords(&records);
  EXPECT_EQ(records.size(), 0);


  traverse.SeekFroward(GenomeLoc("chr1", 0, 10440, 10519));
  traverse.GetRecords(&records);
  EXPECT_EQ(records.size(), 2);
}

TEST(SortingVariantContextWriter, normal) {
  TEST_FILE("vcf/exampleDBSNP.vcf", filename);
  TEST_FILE("vcf/sorting-example.vcf", output_filename);
  VCFWriter vcf_writer(output_filename);
  VCFIndexReader reader_header(filename);
  SortingVariantContextWriter writer(&vcf_writer, 100000);
  VCFHeader* header = writer.GetHeader();
  header->Copy(reader_header.GetHeader());
  writer.WriteHeader();

  auto fun = [&filename, &writer]() {
    VCFTextReader reader(filename);
    while (reader.HasNext()) {
      std::unique_ptr<VariantContext> record(reader.Next());
      writer.Add(record);
    }
  };

  std::vector<std::thread> threads;
  int threads_cnt = 16;
  for (int i = 0; i < threads_cnt; i++) {
    threads.emplace_back(fun);
  }
  for (int i = 0; i < threads_cnt; i++) {
    threads[i].join();
  }
  writer.Close();
}
