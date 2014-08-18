#ifndef QPROFILE_H_
#define QPROFILE_H_

#include<cassert>
#include<fstream>

#include "utils.h"
#include "sampling.hpp"

class QProfile {
public:
  QProfile();
  QProfile& operator=(const QProfile&);
 
  // qual starts from 0, 33 is already deducted
  double getProb(int qual, int ref_base, int read_base) {
    return p[qual][ref_base][read_base];
  }

  void update(int qual, int ref_base, int read_base, double frac) {
    p[qual][ref_base][read_base] += frac;
  }

  void init();
  void collect(const QProfile*);
  void finish();
    
  void read(std::ifstream&);
  void write(std::ofstream&);
 
  char simulate(Sampler* sampler, int qual, int ref_base) {
    return getCharacter(sampler->sample(pc[qual][ref_base], NCODES));
  }
  
  void startSimulation();
  void finishSimulation();
  
private:
  static const int NCODES = 5; // number of possible codes
  static const int SIZE = 100;
  
  double p[SIZE][NCODES][NCODES]; // p[q][r][c] = p(c|r,q)
  
  double (*pc)[NCODES][NCODES]; // for simulation
};

#endif /* QPROFILE_H_ */
