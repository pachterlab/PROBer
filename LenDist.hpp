#ifndef LENDIST_H_
#define LENDIST_H_

#include<cassert>
#include<fstream>
#include<algorithm>

#include "utils.h"
#include "sampling.hpp"

class LenDist {
public:
  /*
    @param minL       minimum length
    @param maxL       maximum length
    @param hasCount   This is a mate length distribution in the preprocessing step, we need to record the counts for different lengths in unalignable reads
   */
  LenDist(int minL = 1, int maxL = 1000, bool hasCount = false);
  ~LenDist();  
  LenDist& operator=(const LenDist&);
  void setAsNormal(double, double);

  int getMinL() const { return lb + 1; }
  int getMaxL() const { return ub; }

  void locate_len_lr(int);

  int get_len_l() const { return len_l; }
  int get_len_r() const { return len_r; }

  double getProb(int len, int upper_bound = MAXV) const {
    if (len <= lb || len > ub || len > upper_bound) return 0.0;
    return (len < upper_bound ? pmf[len - lb] : cdf[span] - cdf[len - 1 - lb]);
  }
  
  // get cumulative probabilities
  double getCProb(int len, int upper_bound = MAXV) const {
    if (len <= lb) return 0.0;
    return ((len < ub && len < upper_bound) ? cdf[len - lb] : 1.0);
  }

  // get partial sum of probability mass
  double getSProb(int from, int to, int upper_bound = MAXV) const {
    if (from > to || to <= lb) return 0.0;
    return ((to < ub && to < upper_bound) ? cdf[to - lb] : 1.0) - cdf[std::max(from - 1 - lb, 0)];
  }

  //len : mate/fragment length
  //refL : reference sequence length, in fact, this is totLen for global length distribution
  double getAdjustedProb(int len, int refL) const {
    if (len <= lb || len > ub || refL <= lb) return 0.0;
    double denom = cdf[std::min(ub, refL) - lb];
    assert(!isLongZero(denom));
    return pmf[len - lb] / denom;
  }
  
  //len : length threshold, any length <= len should be calculated
  //refL : reference sequence length
  double getAdjustedCumulativeProb(int len, int refL) const {
    if (len <= lb || refL <= lb) return 0.0;
    if (len >= ub || len >= refL) return 1.0;
    double denom = cdf[std::min(ub, refL) - lb];
    assert(!isLongZero(denom));
    return cdf[len - lb] / denom;
  }
  
  // We assume the end case will not affect 
  void update(int len, double frac) {
    assert(len > lb && len <= ub);
    pmf[len - lb] += frac;
  }
  
  // Update counts vector
  void updateC(int len) { counts[len] += 1.0; }

  void init();
  void collect(const LenDist* o);  
  void finish();

  // Calculate log probability for unalignable reads
  double calcLogP();

  // Calculate the mean length
  double calcMean() {
    double mean = 0.0;
    for (int i = 1; i <= span; i++) mean += i * pmf[i];
    mean += lb;
    return mean;
  }

  void read(std::ifstream&);
  void write(std::ofstream&);

  int simulate(Sampler* sampler, int upper_bound = MAXV) {
    assert(upper_bound > lb);
    int len = lb + 1 + sampler->sample(cdf + 1, span);
    if (upper_bound < len) len = upper_bound;
    return len;
  }

  int simulateAdjusted(Sampler* sampler, int refL) {
    int dlen;
    if (refL <= lb || cdf[(dlen = std::min(ub, refL) - lb)] <= 0.0) return -1;
    return lb + 1 + sampler->sample(cdf + 1, dlen);
  }

  // Trim zero entries from left and right ends
  void trim();

private:
  static const int MAXV = 2147483647; // the maximum int value

  int lb, ub, span; // (lb, ub]
  double *pmf, *cdf; // probability mass function, cumulative density function

  double *counts; // Counts of different read lengths from unalignable reads

  int len_l, len_r; // [len_l, len_r] gives the lengths used for integration
};

#endif /* LENDIST_H_ */
