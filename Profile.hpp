#ifndef PROFILE_H_
#define PROFILE_H_

#include<fstream>

#include "utils.h"
#include "sampling.hpp"

class Profile {
public:
  Profile(int = 1000);
  ~Profile();
  Profile& operator=(const Profile&);

  double getProb(int pos, int ref_base, int read_base) {
    return p[pos][ref_base][read_base];
  }

  void update(int pos, int ref_base, int read_base, double frac) {
    p[pos][ref_base][read_base] += frac;
  }

  void init();
  void collect(const Profile*);
  void finish();

  void read(std::ifstream&);
  void write(std::ofstream&);

  char simulate(Sampler* sampler, int pos, int ref_base) {
    return getCharacter(sampler->sample(pc[pos][ref_base], NCODES));
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
