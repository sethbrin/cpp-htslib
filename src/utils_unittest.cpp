//
// Created by zp on 12/8/16.
//

#include "easehts/unittest.h"
#include "easehts/utils.h"

#include <gtest/gtest.h>

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
