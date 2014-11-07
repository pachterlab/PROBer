#include<cstring>
#include<cassert>
#include<vector>
#include<string>
#include<fstream>

#include "utils.h"
#include "Markov.hpp"

const std::vector<int> Markov::chr2state = init_chr2state();

std::vector<int> Markov::init_chr2state() {
  std::vector<int> vec(128, -1);
  vec['M'] = vec['='] = vec['X'] = M;
  vec['I'] = I;
  vec['D'] = D;
  return vec;
}

Markov::Markov() {
  // Initialize parameters
  P_start[M] = 0.99; P_start[I] = P_start[D] = 0.005;
  P_trans[M][M] = 0.99; P_trans[M][I] = P_trans[M][D] = 0.005;
  P_trans[I][M] = 0.99; P_trans[I][I] = P_trans[I][D] = 0.005;
  P_trans[D][M] = 0.99; P_trans[D][I] = P_trans[D][D] = 0.005;

  probI[0] = probI[1] = probI[2] = probI[3] = (1.0 - 1e-5) / 4.0;
  probI[4] = 1e-5;

  P_start_sim = NULL;
  P_trans_sim = NULL;
  probI_sim = NULL;
}

void Markov::init() {
  memset(P_start, 0, sizeof(P_start));
  memset(P_trans, 0, sizeof(P_trans));
  memset(probI, 0, sizeof(probI));
}

void Markov::collect(const Markov* o) {
  for (int i = 0; i < NSTATES; ++i) P_start[i] += o->P_start[i];
  for (int i = 0; i < NSTATES; ++i) 
    for (int j = 0; j < NSTATES; ++j) 
      P_trans[i][j] += o->P_trans[i][j];
  for (int i = 0; i < NCODES; ++i) probI[i] += o->probI[i];
}

void Markov::finish() {
  double sum = 0.0;

  for (int i = 0; i < NSTATES; ++i) sum += P_start[i];
  assert(sum > 0.0);
  for (int i = 0; i < NSTATES; ++i) P_start[i] /= sum;

  for (int i = 0; i < NSTATES; ++i) {
    sum = 0.0;
    for (int j = 0; j < NSTATES; ++j) sum += P_trans[i][j];
    if (isZero(sum)) memset(P_trans[i], 0, sizeof(double) * NSTATES);
    else for (int j = 0; j < NSTATES; ++j) P_trans[i][j] /= sum;
  }

  sum = 0.0;
  for (int i = 0; i < NCODES; ++i) sum += probI[i];
  if (isZero(sum)) memset(probI, 0, sizeof(probI));
  else for (int i = 0; i < NCODES; ++i) probI[i] /= sum;
}

void Markov::read(std::ifstream& fin) {
  std::string line;
  while (getline(fin, line)) {
    if (line.substr(0, 7) == "#Markov") break;
  }
  assert(fin.good());
  
  int nstates, ncodes;
  fin>> nstates>> ncodes;
  assert(fin.good() && (nstates == NSTATES) && (ncodes == NCODES));

  for (int i = 0; i < NSTATES; ++i) assert(fin>> P_start[i]);
  
  for (int i = 0; i < NSTATES; ++i) 
    for (int j = 0; j < NSTATES; ++j) assert(fin>> P_trans[i][j]);

  for (int i = 0; i < NCODES; ++i) assert(fin>> probI[i]);
  
  getline(fin, line);
}

void Markov::write(std::ofstream& fout) {
  fout<< "#Markov, format: NSTATES NCODES; P_start; P_trans; probI"<< std::endl;
  fout<< NSTATES<< '\t'<< NCODES<< std::endl;

  for (int i = 0; i < NSTATES - 1; ++i) fout<< P_start[i]<< '\t';
  fout<< P_start[NSTATES - 1]<< std::endl;

  for (int i = 0; i < NSTATES; ++i) {
    for (int j = 0; j < NSTATES - 1; ++j) fout<< P_trans[i][j]<< '\t';
    fout<< P_trans[i][NSTATES - 1]<< std::endl;
  }

  for (int i = 0; i < NCODES - 1; ++i) fout<< probI[i]<< '\t';
  fout<< probI[NCODES - 1]<< std::endl << std::endl;
}

void Markov::startSimulation() {
  P_start_sim = new double[NSTATES];
  memcpy(P_start_sim, P_start, sizeof(P_start));
  for (int i = 1; i < NSTATES; ++i) P_start_sim[i] += P_start_sim[i - 1];

  P_trans_sim = new double[NSTATES][NSTATES];
  memcpy(P_trans_sim, P_trans, sizeof(P_trans));
  for (int i = 0; i < NSTATES; ++i) 
    for (int j = 1; j < NSTATES; ++j) P_trans_sim[i][j] += P_trans_sim[i][j - 1];

  probI_sim = new double[NCODES];
  memcpy(probI_sim, probI, sizeof(probI));
  for (int i = 1; i < NCODES; ++i) probI_sim[i] += probI_sim[i - 1];    
}

void Markov::finishSimulation() {
  delete[] P_start_sim;
  delete[] P_trans_sim;
  delete[] probI_sim;
}
