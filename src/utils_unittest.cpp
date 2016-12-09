//
// Created by zp on 12/8/16.
//

#include "easehts/unittest.h"
#include "easehts/utils.h"

#include <gtest/gtest.h>
#include <vector>

using namespace ncic::easehts::utils;

TEST(Random, NextInt) {
  ThreadLocalRandom r(1111111L);

  EXPECT_EQ(r.NextInt(), -875981731);
  EXPECT_EQ(r.NextInt(), 199965190);
  EXPECT_EQ(r.NextInt(), -306597357);
}

TEST(Random, NextIntWithParameter) {
  ThreadLocalRandom r(1111111L);

  EXPECT_EQ(r.NextInt(4), 3);
  EXPECT_EQ(r.NextInt(5), 0);
  EXPECT_EQ(r.NextInt(1024), 950);
  EXPECT_EQ(r.NextInt(2049), 1996);
}

TEST(SampleIndicesWithoutReplacement, simple) {
  ThreadLocalRandom rand(47382911LL);
  std::vector<int> res = SampleIndicesWithoutReplacemement(9, 5, rand);

  EXPECT_EQ(res[0], 2);
  EXPECT_EQ(res[1], 3);
  EXPECT_EQ(res[2], 4);
  EXPECT_EQ(res[3], 1);
  EXPECT_EQ(res[4], 0);
}
