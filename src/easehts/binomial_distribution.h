//
// Created by zp on 11/25/16.
//

#ifndef EASEHTSLIB_BINOMIAL_DISTRIBUTION_H_
#define EASEHTSLIB_BINOMIAL_DISTRIBUTION_H_

#include <cmath>

namespace ncic {
namespace easehts {

class SaddlePointExpansion {
 private:
  const static double HALF_LOG_2_PI;
  // exact Stirling expansion error for certain values.
  const static double EXACT_STIRLING_ERRORS[32];
  const static double EPS;
 public:
  /**
   * Compute the error of Stirling's series at the given value.
   *
   * @param z the value.
   * @return the Striling's series error.
   */
  static double GetStirlingError(double z) {
    double ret;
    if (z < 15.0) {
      double z2 = 2.0 * z;
      if (std::abs(std::floor(z2) - z2) < SaddlePointExpansion::EPS) {
        ret = SaddlePointExpansion::EXACT_STIRLING_ERRORS[static_cast<int>(z2)];
      } else {
        ret = lgamma(z + 1.0) - (z + 0.5) * std::log(z)
            + z - SaddlePointExpansion::HALF_LOG_2_PI;
      }
    } else {
      double z2 = z * z;
      ret = (0.083333333333333333333 -
             (0.00277777777777777777778 -
              (0.00079365079365079365079365 -
               (0.000595238095238095238095238 -
                0.0008417508417508417508417508 /
                                                    z2) / z2) / z2) / z2) / z;
    }

    return ret;
  }

  /**
   * A part of the deviance portion of the saddle point approximation.
   *
   * @param x  the x value.
   * @param mu the average.
   * @return a part of the deviance.
   */
  static double GetDeviancePart(double x, double mu) {
    double ret;
    if (std::abs(x - mu) < 0.1 * (x + mu)) {
      double d = x - mu;
      double v = d / (x + mu);
      double s1 = v * d;
      double s = NAN;
      double ej = 2.0 * x * v;
      v = v * v;
      int j = 1;
      while (s1 != s) {
        s = s1;
        ej *= v;
        s1 = s + ej / ((j * 2) + 1);
        ++j;
      }
      ret = s1;
    } else {
      ret = x * std::log(x / mu) + mu - x;
    }
    return ret;
  }

  /**
   * Compute the PMF for a binomial distribution using the saddle point
   * expansion.
   *
   * @param x the value at which the probability is evaluated.
   * @param n the number of trials.
   * @param p the probability of success.
   * @param q the probability of failure (1 - p).
   * @return log(p(x)).
   */
  static double LogBinomialProbability(int x, int n, double p, double q) {
    double ret;
    if (x == 0) {
      if (p < 0.1) {
        ret = -GetDeviancePart(n, n * q) - n * p;
      } else {
        ret = n * std::log(q);
      }
    } else if (x == n) {
      if (q < 0.1) {
        ret = -GetDeviancePart(n, n * p) - n * q;
      } else {
        ret = n * std::log(p);
      }
    } else {
      ret = GetStirlingError(n) - GetStirlingError(x) -
          GetStirlingError(n - x) - GetDeviancePart(x, n * p) -
          GetDeviancePart(n - x, n * q);
      double f = (M_PI * 2 * x * (n - x)) / n;
      ret = -0.5 * std::log(f) + ret;
    }
    return ret;
  }
};


class BinomialDistribution
{
 public:
  BinomialDistribution(int trials, double p)
      : number_of_trials_(trials),
      probability_of_success_(p) {}

   /**
    * For this distribution, X, this method returns P(X = x).
    *
    * @param x the value at which the PMF is evaluated.
    * @return PMF for this distribution.
    */
  double Probability(int x) {
    double ret;
    if (x < 0 || x > number_of_trials_) {
      ret = 0.0;
    } else {
      ret = std::exp(SaddlePointExpansion::LogBinomialProbability(x,
          number_of_trials_, probability_of_success_,
          1.0 - probability_of_success_));
    }
    return ret;
  }


  //double cdf(double x);
 private:
  // The number of trials.
  int number_of_trials_;
  // The probability of success.
  double probability_of_success_;
};

} // easehts
} // ncic

#endif
