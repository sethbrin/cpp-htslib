//
// Created by zp on 11/23/16.
//

#ifndef MUTECT_POWER_CALCULATOR_H_
#define MUTECT_POWER_CALCULATOR_H_

#include <easehts/noncopyable.h>
#include <easehts/binomial_distribution.h>

#include <cmath>
#include <functional>
#include <unordered_map>
#include <vector>

namespace ncic {
namespace mutect {

struct PowerCacheKey {
  int n;
  double delta;

  PowerCacheKey(int n_, double delta_) {
    n = n_;
    delta = delta_;
  }

  bool operator==(const PowerCacheKey& rhs) const {
    if (&rhs == this) return true;
    const int EPS = 1e-8;
    if (n != rhs.n) return false;
    if (std::abs(delta - rhs.delta) > EPS) return false;
    return true;
  }
};

} // mutect
} // ncic

// the hash function for PowerCacheKey
namespace std {
template<> struct hash<ncic::mutect::PowerCacheKey>
{
  size_t operator()(ncic::mutect::PowerCacheKey const& rhs) const
  {
    size_t const h1(std::hash<int>{}(rhs.n));
    size_t const h2(std::hash<double>{}(rhs.delta));
    return h1 ^ h2;
  }
};
}

namespace ncic {
namespace mutect {

class AbstractPowerCalculator {
 public:
  explicit AbstractPowerCalculator() {}
  AbstractPowerCalculator(double constant_eps, double constant_lod_threshold)
    : constant_eps_(constant_eps),
    constant_lod_threshold_(constant_lod_threshold) {}

  static double CalculateLogLikelihood(int depth, int alts, double eps, double f) {
    double a = (depth - alts) * std::log10(f * eps + (1 - f) * (1 - eps));
    double b = (alts) * std::log10(f * (1 - eps) + (1 - f) * eps);

    return (a + b);
  }

 protected:
  std::unordered_map<PowerCacheKey, double> cache_;
  double constant_eps_;
  double constant_lod_threshold_;
};

class NormalPowerCalculator : public AbstractPowerCalculator {
 public:
  explicit NormalPowerCalculator() {}
  NormalPowerCalculator(double constant_eps, double constant_lod_threshold)
    : NormalPowerCalculator(constant_eps, constant_lod_threshold, true) {}

  NormalPowerCalculator(double constant_eps, double constant_lod_threshold, bool enable_smoothing)
    : AbstractPowerCalculator(constant_eps, constant_lod_threshold),
    constant_enable_smoothing_(enable_smoothing) {}

  double CachingPowerCalculation(int n) {
    PowerCacheKey key(n, 0.5);
    if (cache_.find(key) == cache_.end()) {
      double power = CalculatePower(n, constant_eps_,
                                    constant_lod_threshold_, constant_enable_smoothing_);
      cache_[key] = power;
      return power;
    }
    return cache_[key];
  }

  static double CalculateNormalLod(int depth, int alts, double eps) {
    return (CalculateLogLikelihood(depth, alts, eps, 0) -
            CalculateLogLikelihood(depth, alts, eps, 0.5));
  }

  static double CalculatePower(int depth, double eps, double lod_threshold, bool enable_smoothing) {
    if (depth == 0) return 0;

    // Calculate the probalility of each configuration
    // NOTE: reorder from tumor case to be # of ref instead of # of alts
    easehts::BinomialDistribution binom(depth, eps);
    std::vector<double> p(depth + 1);
    std::vector<double> lod(depth + 1);
    int k = -1;
    for (int i = 0; i < p.size(); i++) {
      p[i] = binom.Probability(depth - i);
      lod[i] = CalculateNormalLod(depth, depth - i, eps);

      if (lod[i] >= lod_threshold && k == -1) {
        k = i;
      }
    }

    // if no depth meets the lod score, the power is zero
    if (k == -1) return 0;

    double power = 0;
    // here we correct for the fact that the exact lod threshold is likely somewhere between
    // the k and k-1 bin, so we prorate the power from that bin
    // the k and k-1 bin, so we prorate the power from that bin
    // if k==0, it must be that lodThreshold == lod[k] so we don't have to make this correction
    if (enable_smoothing && k > 0) {
      double x = 1 - (lod_threshold - lod[k-1]) / (lod[k] - lod[k-1]);
      power = x * p[k-1];
    }
    for (int i = k; i < p.size(); i++) {
      power += p[i];
    }

    return power;

  }

 private:
  bool constant_enable_smoothing_;
};

class TumorPowerCalculator : public AbstractPowerCalculator {
 public:
  explicit TumorPowerCalculator() {}
  TumorPowerCalculator(double constant_eps, double constant_lod_threshold,
                       double constant_contamination)
    : TumorPowerCalculator(constant_eps, constant_lod_threshold,
                            constant_contamination, true) {}

  TumorPowerCalculator(double constant_eps, double constant_lod_threshold,
                       double constant_contamination, bool enable_smoothing)
    : AbstractPowerCalculator(constant_eps, constant_lod_threshold),
    constant_enable_smoothing_(enable_smoothing),
    constant_contamination_(constant_contamination){}

  double CachingPowerCalculation(int n, double delta) {
    PowerCacheKey key(n, delta);
    if (cache_.find(key) == cache_.end()) {
      double power = CalculatePower(n, constant_eps_,
                                    constant_lod_threshold_, delta,
                                    constant_contamination_,
                                    constant_enable_smoothing_);
      cache_[key] = power;
      return power;
    }
    return cache_[key];
  }

  static double CalculateNormalLod(int depth, int alts, double eps, double contam) {
    double f = static_cast<double>(alts) / static_cast<double>(depth);
    return (CalculateLogLikelihood(depth, alts, eps, f) -
            CalculateLogLikelihood(depth, alts, eps, std::min(f, contam)));
  }

  static double CalculatePower(int depth, double eps, double lod_threshold,
                               double delta, double contam, bool enable_smoothing) {
    if (depth == 0) return 0;

    // Calculate the probalility of each configuration
    // NOTE: reorder from tumor case to be # of ref instead of # of alts
    double p_alt_given_e_delta = delta * (1 - eps) + (1 - delta) * eps;
    easehts::BinomialDistribution binom(depth, p_alt_given_e_delta);
    std::vector<double> p(depth + 1);
    std::vector<double> lod(depth + 1);
    int k = -1;
    for (int i = 0; i < p.size(); i++) {
      p[i] = binom.Probability(i);
      lod[i] = CalculateNormalLod(depth, i, eps, contam);

      if (lod[i] >= lod_threshold && k == -1) {
        k = i;
      }
    }

    // if no depth meets the lod score, the power is zero
    if (k == -1) return 0;

    double power = 0;
    // here we correct for the fact that the exact lod threshold is likely somewhere between
    // the k and k-1 bin, so we prorate the power from that bin
    // the k and k-1 bin, so we prorate the power from that bin
    // if k==0, it must be that lodThreshold == lod[k] so we don't have to make this correction
    if (enable_smoothing && k > 0) {
      double x = 1 - (lod_threshold - lod[k-1]) / (lod[k] - lod[k-1]);
      power = x * p[k-1];
    }
    for (int i = k; i < p.size(); i++) {
      power += p[i];
    }

    return power;

  }

 private:
  bool constant_enable_smoothing_;
  double constant_contamination_;
};


} // mutect
} // ncic

#endif
