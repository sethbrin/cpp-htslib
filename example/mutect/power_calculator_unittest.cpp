#include "power_calculator.h"

#include <easehts/unittest.h>
#include <gtest/gtest.h>

#include <unordered_map>
#include <cmath>

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

double Q30_EPS = std::pow(10, (-3)) / 3;
double LOD = 6.3;
double wiggle = 0.0001;

TEST(AbstractPowerCalculator, CalculateLogLikelihood) {

  EXPECT_NEAR(NormalPowerCalculator::CalculateLogLikelihood(50, 5, 0.001, 0.1), -7.059164, wiggle);
  EXPECT_NEAR(NormalPowerCalculator::CalculateLogLikelihood(40, 5, 0.001, 0.1), -6.597727, wiggle);
  EXPECT_NEAR(NormalPowerCalculator::CalculateLogLikelihood(40, 10, 0.001, 0.1), -11.34971, wiggle);
  EXPECT_NEAR(NormalPowerCalculator::CalculateLogLikelihood(40, 10, 0.01, 0.1), -11.15482, wiggle);
  EXPECT_NEAR(NormalPowerCalculator::CalculateLogLikelihood(40, 10, 0.01, 0.2), -9.866713, wiggle);
}

TEST(TumorPowerCalculator, CalculatePower) {
  EXPECT_NEAR(TumorPowerCalculator::CalculatePower(50, Q30_EPS, LOD, 0.15, 0, true),
              0.975342, wiggle);
  EXPECT_NEAR(TumorPowerCalculator::CalculatePower(20, Q30_EPS, LOD, 0.15, 0, true), 0.6364833, wiggle);
  EXPECT_NEAR(TumorPowerCalculator::CalculatePower(10, Q30_EPS, LOD, 0.15, 0, true), 0.3167154, wiggle);
  EXPECT_NEAR(TumorPowerCalculator::CalculatePower(50, Q30_EPS, LOD, 0.35, 0, true), 0.9999994, wiggle);
  EXPECT_NEAR(TumorPowerCalculator::CalculatePower(50, Q30_EPS, LOD, 0.05, 0, true), 0.3893047, wiggle);
  EXPECT_NEAR(TumorPowerCalculator::CalculatePower(50, Q30_EPS, 4.3, 0.05, 0, true), 0.6068266, wiggle);
  EXPECT_NEAR(TumorPowerCalculator::CalculatePower(50, Q30_EPS, 4.3, 0.05, 0.02, true), 0.0039135, wiggle);
}

TEST(TumorPowerCalculator, CachingPowerCalculation) {
  TumorPowerCalculator pc(Q30_EPS, LOD, 0);
  EXPECT_NEAR(pc.CachingPowerCalculation(50, 0.15), 0.975342, wiggle);
  EXPECT_NEAR(pc.CachingPowerCalculation(50, 0.15), 0.975342, wiggle);
  EXPECT_NEAR(pc.CachingPowerCalculation(50, 0.15), 0.975342, wiggle);
}
