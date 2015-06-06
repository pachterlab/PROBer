#ifndef NOISEPROFILE_H_
#define NOISEPROFILE_H_

#include<cstring>
#include<string>
#include<fstream>

#include "utils.h"
#include "sampling.hpp"
#include "SEQstring.hpp"

class NoiseProfile {
public:
  NoiseProfile(bool hasCount = false);
  ~NoiseProfile();
  
  void updateC(const SEQstring& seq) {
    int len = seq.getLen();
    for (int i = 0; i < len; ++i) {
      ++c[seq.baseCodeAt(i)];
    }
  }

  void writeC(std::ofstream& fout);
  void readC(std::ifstream& fin);

  void calcInitParams();

  double getProb(const SEQstring& seq) {
    double prob = 1.0;
    int len = seq.getLen();
    
    for (int i = 0; i < len; ++i) {
      prob *= p[seq.baseCodeAt(i)];
    }
    
    return prob;
  }
  
  void update(const SEQstring& seq, double frac) {
    int len = seq.getLen();
    for (int i = 0; i < len; ++i) {
      p[seq.baseCodeAt(i)] += frac;
    }
  }

  void init();
  void collect(const NoiseProfile* o);
  void finish();
  
  double calcLogP();

  void read(std::ifstream& fin);
  void write(std::ofstream& fout);

  void simulate(Sampler* sampler, int len, std::string& readseq) {
    readseq.assign(len, 0);
    for (int i = 0; i < len; ++i) {
      readseq[i] = code2base[sampler->sample(pc, NCODES)];
    }
  }
  
  void startSimulation();
  void finishSimulation();
  
private:
  static const int NCODES = 5;

  double p[NCODES];
  double *c; // counts in N0

  double *pc; // for simulation
};

#endif /* NOISEPROFILE_H_ */
