#ifndef ORIENTATION_H_
#define ORIENTATION_H_

#include<cstring>
#include<cassert>
#include<string>
#include<fstream>

#include "sampling.hpp"

class Orientation {
 public:
  Orientation(double probF = 0.5) {
    prob[0] = probF;
    prob[1] = 1.0 - probF;

    eprob[0] = eprob[1] = 0.0;
  }
  
  Orientation& operator= (const Orientation& rv) {
    if (this == &rv) return *this;
    memcpy(prob, rv.prob, sizeof(rv.prob));

    memcpy(eprob, rv.eprob, sizeof(rv.eprob));

    return *this;
  }
  
  //dir : +/-
  double getProb(char dir) { return prob[dir == '+' ? 0 : 1]; }
  
  void init() { eprob[0] = eprob[1] = 0.0; }

  //dir must be either + or -
  void update(char dir) { ++eprob[dir == '+' ? 0 : 1]; }

  void collect(const Orientation* o) {
    eprob[0] += o->eprob[0];
    eprob[1] += o->eprob[1];
  }

  void finish() { 
    double sum = eprob[0] + eprob[1];
    eprob[0] /= sum;
    eprob[1] /= sum;
  }

  void read(std::ifstream& fin) {
    std::string line;
    while (getline(fin, line)) {
      if (line.substr(0, 12) == "#Orientation") break;
    }
    assert(fin.good()); // even eofbit is not set

    fin>> prob[0]>> eprob[0]; 
    assert(fin); // test if read is successful, it's OK to have the eofbit set
    getline(fin, line);

    prob[1] = 1.0 - prob[0];
    eprob[1] = 1.0 - eprob[0];
  }
  
  void write(std::ofstream& fout) {
    fout<< "#Orientation, the set forward probability and estimated forward probability"<< std::endl;
    fout.precision(10);
    fout.setf(0, std::ios::floatfield);
    fout<< prob[0]<< '\t'<< eprob[0]<< std::endl << std::endl;
  }
  
  char simulate(Sampler* sampler) { return (sampler->random() < prob[0] ? '+' : '-'); }
  
 private:  
  double prob[2]; // 0 + 1 -
  double eprob[2]; // estimated probabilities
};

#endif /* ORIENTATION_H_ */
