#ifndef MARKOV_H_
#define MARKOV_H_

/**
 * This class models a Markov chain for insertions/deletions. Because insertions and deletions are rare for Illumina reads, the insertions/deletions are assumed to be due to an incomplete 
 * knowledge of the reference. 
 */

#include<vector>
#include<fstream>

#include "utils.h"
#include "sampling.hpp"

class Markov {
public:
  Markov();
  Markov& operator=(const Markov&);

  double getProb(char a) { return P_start[chr2state[a]]; }
  double getProb(char a, char b) { return P_trans[chr2state[a]][chr2state[b]]; }
  
  double getIBaseProb(int code) { return probI[code]; }

  void update(char a, double frac) { P_start[chr2state[a]] += frac; }
  void update(char a, char b, double frac) { P_trans[chr2state[a]][chr2state[b]] += frac; }

  void updateIBase(int code, double frac) { probI[code] += frac; }

  void init();
  void collect(const Markov*);
  void finish();

  void read(std::ifstream&);
  void write(std::ofstream&);

  char simulate(Sampler *sampler) {
    switch(sampler->sample(P_start_sim, NSTATES)) {
    case M : return 'M';
    case I : return 'I';
    case D : return 'D';
    default : assert(false);
    }
  }

  char simulate(Sampler *sampler, char a) {
    switch(sampler->sample(P_trans_sim[chr2state[a]], NSTATES)) {
    case M : return 'M';
    case I : return 'I';
    case D : return 'D';
    default : assert(false);
    }
  }

  char simulateIBase(Sampler *sampler) {
    return code2base[sampler->sample(probI_sim, NCODES)];
  }

  void startSimulation();
  void finishSimulation();

private:
  static const int NSTATES = 3;
  static const int M = 0;
  static const int I = 1;
  static const int D = 2;

  static const std::vector<int> chr2state;
  static std::vector<int> init_chr2state();

  double P_start[NSTATES]; // model the probability of a read starts from a state
  double P_trans[NSTATES][NSTATES];

  static const int NCODES = 5; // ACGTN
  
  double probI[NCODES]; // The probability of generating a base given the state is I

  // for simulation
  double *P_start_sim;
  double (*P_trans_sim)[NSTATES];
  double *probI_sim;
};

#endif 
