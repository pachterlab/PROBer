#include<cmath>
#include<cassert>
#include<string>
#include<fstream>

#include "MateLenDist.hpp"

MateLenDist::MateLenDist() {
  lb = ub = span = 0;
  logp = 0.0;
  
  pmf.clear();
  cdf.clear();
  noise_counts.clear();
}

void MateLenDist::finish() {
  double sum = 0.0;

  assert(ub > 0);
  for (lb = 1; pmf[lb] == 0.0; ++lb);
  for (int i = lb; i <= ub; ++i) {
    sum += pmf[i];
    pmf[i - lb] = pmf[i];
  }
  span = ub - lb + 1;
  pmf.resize(span);
  for (int i = 0; i < span; ++i) pmf[i] /= sum;

  cdf.resize(span);
  for (int i = 0; i < span; ++i) cdf[i] = pmf[i] + (i > 0 ? cdf[i - 1] : 0.0);

  // Calculate logp
  logp = 0.0;
  for (int i = 0; i < span; ++i) 
    if (noise_counts[i + lb] > 0.0) logp += noise_counts[i + lb] * log(pmf[i]);
}

void MateLenDist::read(std::ifstream& fin) {
  std::string line;
  while (getline(fin, line)) {
    if (line.substr(0, 12) == "#MateLenDist") break;
  }
  assert(fin.good());

  assert(fin>> lb>> ub>> span);
  pmf.resize(span, 0.0);
  cdf.resize(span, 0.0);
  for (int i = 0; i < span; ++i) {
    assert(fin>> pmf[i]);
    cdf[i] = pmf[i];
    if (i > 0) cdf[i] += cdf[i - 1];
  }
}

void MateLenDist::write(std::ofstream& fout) {
  fout<< "#MateLenDist, format: lb ub span; [lb, ub], span = ub - lb + 1, probability mass function values"<< std::endl;
  fout<< lb<< ub<< span<< std::endl;  
  for (int i = 0; i < span - 1; ++i) fout<< pmf[i]<< '\t';
  fout<< pmf[span - 1]<< std::endl<< std::endl;
}
