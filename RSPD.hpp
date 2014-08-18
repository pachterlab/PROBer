#ifndef RSPD_H_
#define RSPD_H_

#include<cassert>
#include<fstream>

#include "sampling.hpp"

class RSPD {
public:
  RSPD(bool estRSPD = false) {
    this->estRSPD = estRSPD;
    assert(!estRSPD); 
  }
  
  RSPD& operator=(const RSPD& rv) {
    if (this == &rv) return *this;
    estRSPD = rv.estRSPD;
    return *this;
  }

  bool estimateRSPD() const { return estRSPD; }

  double getProb(int fpos, int effL) { 
    assert(!estRSPD);
    return 1.0 / effL;
  }

  void init() { assert(false); }
  void collect(RSPD* o) { assert(false); }
  void finish() { assert(false); }

  void read(std::ifstream&);
  void write(std::ofstream&);

  int simulate(Sampler* sampler, int effL) {
    return int(sampler->random() * effL);
  }
  
private:
  bool estRSPD;
};

#endif /* RSPD_H_ */
