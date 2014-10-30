#ifndef PROFILE_H_
#define PROFILE_H_

#include<fstream>

#include "utils.h"
#include "sampling.hpp"

class Profile {
public:
  Profile(int maxL = 1000);
  ~Profile();

  double getProb(int pos, int ref_base, int read_base) {
    return p[pos][ref_base][read_base];
  }

  void update(int pos, int ref_base, int read_base, double frac) {
    p[pos][ref_base][read_base] += frac;
  }

  void init();
  void collect(const Profile* o);
  void finish();

  void read(std::ifstream& fin);
  void write(std::ofstream& fout);

  char simulate(Sampler* sampler, int pos, int ref_base) {
    return code2base[sampler->sample(pc[pos][ref_base], NCODES)];
  }

  void startSimulation();
  void finishSimulation();
  
private:
  static const int NCODES = 5;
  
  int proLen; // profile length
  int size; // # of items in p;
  double (*p)[NCODES][NCODES]; //profile matrices
  
  double (*pc)[NCODES][NCODES]; // for simulation
};

#endif /* PROFILE_H_ */
