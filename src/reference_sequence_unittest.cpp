//
// Created by zp on 11/18/16.
//

#include "easehts.h"
#include "reference_sequence.h"

#include "unittest.h"

#include <string.h>
#include <gtest/gtest.h>

using namespace ncic::easehts;

TEST(FindSequenceDictionaryFileName, AbstractFastaSequenceFile) {
  TEST_FILE("references/auxf.fa", filename);
  EXPECT_STREQ(AbstractFastaSequenceFile::FindSequenceDictionaryFileName(filename).c_str(), "");

  TEST_FILE("references/references/human_g1k_v37.fasta", b37);
  TEST_FILE("references/references/human_g1k_v37", b37expect);
  strcat(b37expect, ".dict");
  EXPECT_STREQ(AbstractFastaSequenceFile::FindSequenceDictionaryFileName(b37).c_str(), b37expect);
}

TEST(IsIndexed, IndexedFastaSequenceFile) {
  TEST_FILE("references/references/human_g1k_v37.fasta", b37);
  IndexedFastaSequenceFile reference(b37);

  EXPECT_EQ(reference.IsIndexed(), true);
}

TEST(GetSequenceAt, IndexedFastaSequenceFile) {
  TEST_FILE("references/references/human_g1k_v37.fasta", b37);
  IndexedFastaSequenceFile reference(b37);
  EXPECT_STREQ(reference.GetSequenceAt("1", 0, 1).base, "NN");
}
