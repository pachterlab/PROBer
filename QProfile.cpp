#include<cstring>
#include<cassert>
#include<string>
#include<fstream>

#include "utils.h"
#include "QProfile.hpp"

QProfile::QProfile() {
  memset(p, 0, sizeof(p));
  
  //make initialized parameters
  //ASSUME order of A, C, G, T, N
  int N = NCODES - 1;
  double probN = 1e-5;
  double probC, probO; // current, other

  // Probabillity of N is fixed as 1e-5, 
  // probability of generating a base correct is (1 - 1e-5) * (1 - probO), 
  // where probO is converted back from the quality score
  for (int i = 0; i < SIZE; i++) {
    for (int j = 0; j < NCODES - 1; j++) {
      p[i][j][N] = probN;
      
      probO = exp(-i / 10.0 * log(10.0));
      probC = 1.0 - probO;
      probO /= (NCODES - 2);
      
      probC *= (1.0 - probN);
      probO *= (1.0 - probN);
      
      assert(probC >= 0.0 && probO >= 0.0);
      
      for (int k = 0; k < NCODES - 1; k++) {
	if (j == k) p[i][j][k] = probC;
	else p[i][j][k] = probO;
      }
    }
    p[i][N][N] = probN;
    for (int k = 0; k < NCODES - 1; k++)
      p[i][N][k] = (1.0 - probN) / (NCODES - 1);
  }
}

QProfile& QProfile::operator=(const QProfile& rv) {
  if (this == &rv) return *this;
  memcpy(p, rv.p, sizeof(rv.p));
  return *this;
}

void QProfile::init() {
  memset(p, 0, sizeof(p));
}

void QProfile::collect(const QProfile* o) {
  for (int i = 0; i < SIZE; i++)
    for (int j = 0; j < NCODES; j++)
      for (int k = 0; k < NCODES; k++)
	p[i][j][k] += o->p[i][j][k];
}

void QProfile::finish() {
  double sum;
  
  for (int i = 0; i < SIZE; i++) {
    for (int j = 0; j < NCODES; j++) {
      sum = 0.0;
      for (int k = 0; k < NCODES; k++) sum += p[i][j][k];
      if (isZero(sum)) memset(p[i][j], 0, sizeof(double) * NCODES);
      else for (int k = 0; k < NCODES; k++) p[i][j][k] /= sum;
    }
  }
}

void QProfile::read(std::ifstream& fin) {
  std::string line;
  while (getline(fin, line)) {
    if (line.substr(0, 9) == "#QProfile") break;
  }
  assert(fin.good());

  int tmp_size, tmp_ncodes;
  fin>> tmp_size>> tmp_ncodes;
  assert(fin.good() && (tmp_size == SIZE) && (tmp_ncodes == NCODES));

  for (int i = 0; i < SIZE; i++)
    for (int j = 0; j < NCODES; j++)
      for (int k = 0; k < NCODES; k++)
	assert(fin>> p[i][j][k]);

  getline(fin, line);
}

void QProfile::write(std::ofstream& fout) {
  fout<< "#QProfile, format: SIZE NCODES; P[QUAL][REF_BASE][OBSERVED_BASE], SIZE blocks separated by a blank line, each block contains NCODES lines"<< std::endl;
  fout<< SIZE<< '\t'<< NCODES<< std::endl;

  fout.precision(10);
  fout.setf(0, std::ios::floatfield);
  
  for (int i = 0; i < SIZE; i++) {
    for (int j = 0; j < NCODES; j++) {
      for (int k = 0; k < NCODES - 1; k++)
	fout<< p[i][j][k];
      fout<< p[i][j][NCODES - 1]<< std::endl;
    }
    fout<< std::endl;
  }
}

void QProfile::startSimulation() {
  pc = new double[SIZE][NCODES][NCODES];

  for (int i = 0; i < SIZE; i++) {
    for (int j = 0; j < NCODES; j++)
      for (int k = 0; k < NCODES; k++) {
	pc[i][j][k] = p[i][j][k];
	if (k > 0) pc[i][j][k] += pc[i][j][k - 1];
      }
    
    // In case one (Qual, Ref_base) combination is not seen, sharing information from other combinations
    double cp_sum, cp_d, cp_n;
    double p_d, p_o, p_n;
    
    cp_sum = cp_d = cp_n = 0.0;
    for (int j = 0; j < NCODES - 1; j++) {
      cp_sum += pc[i][j][NCODES - 1];
      cp_d += p[i][j][j];
      cp_n += p[i][j][NCODES - 1];
    }
    
    if (isZero(cp_sum)) {
      p_n = 1e-5;
      p_d = (1.0 - exp(-i / 10.0 * log(10.0))) * (1.0 - p_n);
      p_o = (1.0 - p_d - p_n) / (NCODES - 2);
    }
    else {
      p_d = cp_d / cp_sum;
      p_n = cp_n / cp_sum;
      p_o = (1.0 - p_d - p_n) / (NCODES - 2);
    }
    
    // Check if (Qual, j) has no probability for j != N
    for (int j = 0; j < NCODES - 1; j++) {
      if (!isZero(pc[i][j][NCODES - 1])) continue;
      
      for (int k = 0; k < NCODES; k++) {
	if (k == j) pc[i][j][k] = p_d;
	else if (k == NCODES - 1) pc[i][j][k] = p_n;
	else pc[i][j][k] = p_o;
	if (k > 0) pc[i][j][k] += pc[i][j][k - 1];
      }
    }
    
    // Check if (Qual, N) has no probability
    if (isZero(pc[i][NCODES - 1][NCODES - 1])) {
      p_o = (1.0 - p_n) / (NCODES - 1);
      for (int k = 0; k < NCODES; k++) {
	pc[i][NCODES - 1][k] = (k < NCODES - 1 ? p_o : p_n);
	if (k > 0) pc[i][NCODES - 1][k] += pc[i][NCODES - 1][k - 1];
      }
    }    
  }
}

void QProfile::finishSimulation() {
  delete[] pc;
}
