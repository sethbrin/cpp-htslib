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

TEST(RoundNearest, normal) {
  double val = 1.2345678;

  EXPECT_DOUBLE_EQ(1.234568, RoundNearest(val, 6));
  EXPECT_DOUBLE_EQ(1.23457, RoundNearest(val, 5));
  EXPECT_DOUBLE_EQ(1.2346, RoundNearest(val,4));
  EXPECT_DOUBLE_EQ(1.235, RoundNearest(val,3));
  EXPECT_DOUBLE_EQ(1.23, RoundNearest(val,2));
  EXPECT_DOUBLE_EQ(1.2, RoundNearest(val,1));
  EXPECT_DOUBLE_EQ(1, RoundNearest(val,0));
}

TEST(RoundNearestFormat, normal) {
  EXPECT_STREQ(RoundNearestFormat(1.2345678, 6).c_str(), "1.234568");
  EXPECT_STREQ(RoundNearestFormat(1.23456, 6).c_str(), "1.23456");
  EXPECT_STREQ(RoundNearestFormat(1.2345, 6).c_str(), "1.2345");
  EXPECT_STREQ(RoundNearestFormat(1.234, 6).c_str(), "1.234");
  EXPECT_STREQ(RoundNearestFormat(1.23, 6).c_str(), "1.23");
  EXPECT_STREQ(RoundNearestFormat(1.2, 6).c_str(), "1.2");
  EXPECT_STREQ(RoundNearestFormat(1.0, 6).c_str(), "1");
  EXPECT_STREQ(RoundNearestFormat(-236.475073, 6).c_str(), "-236.475073");
  EXPECT_STREQ(RoundNearestFormat(4.116630, 6).c_str(), "4.11663");
  EXPECT_STREQ(RoundNearestFormat(0.28536039716838274, 6).c_str(), "0.28536");
}
