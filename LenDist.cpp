#include<cmath>
#include<cstring>
#include<cassert>
#include<string>
#include<fstream>

#include "boost/math/distributions/normal.hpp"

#include "utils.h"
#include "my_assert.h"
#include "LenDist.hpp"

LenDist::LenDist(int minL, int maxL, bool hasCount) {
  lb = minL - 1;
  ub = maxL;
  span = ub - lb;
  assert(span > 0);

  pmf = cdf = counts = NULL;

  pmf = new double[span + 1];
  cdf = new double[span + 1];
 
  if (hasCount) {
    counts = new double[span + 1];
    memset(pmf, 0, sizeof(double) * (span + 1));
    memset(cdf, 0, sizeof(double) * (span + 1));
    memset(counts, 0, sizeof(double) * (span + 1));
  }
  else {
    //set initial parameters
    pmf[0] = cdf[0] = 0.0;
    for (int i = 1; i <= span; i++) {
      pmf[i] = 1.0 / span;
      cdf[i] = i * 1.0 / span;
    }
  }
}

LenDist::~LenDist() {
  delete[] pmf;
  delete[] cdf;
  if (counts != NULL) delete[] counts;
}

LenDist& LenDist::operator=(const LenDist& rv) {
  if (this == &rv) return *this;
  if (span != rv.span) {
    delete[] pmf;
    delete[] cdf;
    if (counts != NULL) delete[] counts;
    pmf = new double[rv.span + 1];
    cdf = new double[rv.span + 1];
    if (rv.counts != NULL) counts = new double[rv.span + 1];
  }
  lb = rv.lb; ub = rv.ub; span = rv.span;
  memcpy(pmf, rv.pmf, sizeof(double) * (span + 1));
  memcpy(cdf, rv.cdf, sizeof(double) * (span + 1));
  if (rv.counts != NULL) memcpy(counts, rv.counts, sizeof(double) * (span + 1));
  
  return *this;
}

void LenDist::setAsNormal(double mean, double sd) {
  if (isZero(sd)) {
    // Assume all fragment lengths are equal, which is the mean
    ub = int(mean + 0.5); lb = ub - 1; span = 1;
    pmf[0] = cdf[0] = 0.0;
    pmf[1] = cdf[1] = 1.0;
    return;
  }
    
  boost::math::normal norm(mean, sd);
  double sum, old_value, value;

  sum = 0.0;
  old_value = boost::math::cdf(norm, lb + 0.5);
  for (int i = 1; i <= span; i++) {
    value = boost::math::cdf(norm, lb + i + 0.5);
    pmf[i] = value - old_value;
    sum += pmf[i];
    old_value = value;
  }
  assert(!isZero(sum));

  pmf[0] = cdf[0] = 0.0;
  for (int i = 1; i <= span; i++) {
    pmf[i] /= sum;
    cdf[i] = cdf[i - 1] + pmf[i];
  }
  
  trim();
}

// range_approx: the range used for approximation
void LenDist::locate_len_lr(int range_approx) {
  if (span <= range_approx) { len_l = lb + 1; len_r = ub; return; }
  double best = 0.0;
  for (int i = range_approx; i <= span; ++i) 
    if (best < cdf[i] - cdf[i - range_approx]) {
      best = cdf[i] - cdf[i - range_approx];
      len_r = lb + i; len_l = len_r - range_approx + 1;
    }
}

void LenDist::init() {
  memset(pmf, 0, sizeof(double) * (span + 1));
  memset(cdf, 0, sizeof(double) * (span + 1));
}

void LenDist::collect(const LenDist* o) {
  if (lb != o->lb || ub != o->ub) {
    delete[] pmf;
    delete[] cdf;
    lb = o->lb; ub = o->ub; span = o->span;
    pmf = new double[span + 1];
    cdf = new double[span + 1];
    memset(pmf, 0, sizeof(double) * (span + 1));
    memset(cdf, 0, sizeof(double) * (span + 1));
  }
  for (int i = 1; i <= span; i++) {
    pmf[i] += o->pmf[i];
  }
}
 
void LenDist::finish() {
  double sum = 0.0;
  
  for (int i = 1; i <= span; i++) sum += pmf[i];
  general_assert(!isZero(sum), "No valid read to estimate the length distribution!");
  
  for (int i = 1; i <= span; i++) {
    pmf[i] /= sum;
    cdf[i] = cdf[i - 1] + pmf[i];
  }
}

double LenDist::calcLogP() {
  double logp = 0.0;

  assert(counts != NULL);
  for (int i = 1; i <= span; i++)
    if (counts[i] > 0.0) logp += counts[i] * log(pmf[i]);

  return logp;
}

void LenDist::read(std::ifstream& fin) {
  std::string line;
  while (getline(fin, line)) {
    if (line.substr(0, 8) == "#LenDist") break;
  }
  assert(fin.good());

  //release default space first
  delete[] pmf;
  delete[] cdf;

  assert(fin>> lb>> ub>> span);
  pmf = new double[span + 1];
  cdf = new double[span + 1];
  pmf[0] = cdf[0] = 0.0;
  for (int i = 1; i <= span; i++) {
    assert(fin>> pmf[i]);
    cdf[i] = cdf[i - 1] + pmf[i];
  }
}

void LenDist::write(std::ofstream& fout) {
  trim(); // for a succinct representation

  fout<< "#LenDist, format: lb ub span; probability mass function values"<< std::endl;
  fout<< lb<< ub<< span<< std::endl;
  
  fout.precision(10);
  fout.setf(0, std::ios::floatfield);
  for (int i = 1; i < span; i++) fout<< pmf[i];
  fout<< pmf[span]<< std::endl<< std::endl;
}

void LenDist::trim() {
  int newlb, newub;
  
  for (newlb = 1; newlb <= span && isLongZero(pmf[newlb]); ++newlb);
  --newlb;
  for (newub = span; newub > newlb && isLongZero(pmf[newub]); --newub);
  assert(newlb < newub);
  if (newlb == 0 && newub == span) return;

  span = newub - newlb;
  lb += newlb;
  ub = lb + span;

  if (newlb > 0) memmove(pmf + 1, pmf + newlb + 1, sizeof(double) * span);
  double sum = 0.0;
  for (int i = 1; i <= span; ++i) sum += pmf[i];
  for (int i = 1; i <= span; ++i) {
    pmf[i] /= sum;
    cdf[i] = cdf[i - 1] + pmf[i];
  }
}
