#include "power_calculator.h"

#include <easehts/unittest.h>
#include <gtest/gtest.h>

#include <unordered_map>

using namespace ncic::mutect;
TEST(PowerCacheKey, hash) {
  PowerCacheKey key1(1, 2.0);
  PowerCacheKey key2(2, 3.0);
  PowerCacheKey key3(2, 2.0);
  PowerCacheKey key4(2, 3.0);

  std::unordered_map<PowerCacheKey, int> map;
  map[key1] = 2;
  map[key2] = 3;

  EXPECT_TRUE(map.find(key3) == map.end());
  EXPECT_TRUE(map.find(key4) != map.end());
  EXPECT_EQ(map[key4], 3);
}
