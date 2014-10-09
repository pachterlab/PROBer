#ifndef QUALDIST_H_
#define QUALDIST_H_

#include<cstring>
#include<cassert>
#include<string>
#include<fstream>

#include "sampling.hpp"
#include "QUALstring.hpp"

//from 33 to 126 to encode 0 to 93
class QualDist {
public:
  QualDist();
  
  void update(const QUALstring& qual) {
    int len = qual.getLen();

    ++p_init[qual.qualAt(0)];
    for (int i = 1; i < len; ++i) {
      ++p_tran[qual.qualAt(i - 1)][qual.qualAt(i)];
    }
  }

  double getProb(const QUALstring& qual) {
    int len = qual.getLen();
    double prob = 1.0;
    
    prob *= p_init[qual.qualAt(0)];
    for (int i = 1; i < len; ++i) {
      prob *= p_tran[qual.qualAt(i - 1)][qual.qualAt(i)];
    }
    
    return prob;
  }

  double finish();
  
  void read(std::ifstream&);
  void write(std::ofstream&);
  
  void simulate(Sampler* sampler, int len, std::string& qual) {
    int qval, old_qval;

    qual.assign(len, 0);
    qval = sampler->sample(qc_init, SIZE);
    qual[0] = char(qval + 33);
    for (int i = 1; i < len; ++i) {
      old_qval = qval;
      qval = sampler->sample(qc_trans[old_qval], SIZE);
      qual[i] = char(qval + 33);
    }
  }
  
  void startSimulation();
  void finishSimulation();
  
private:
  static const int SIZE = 100;
  
  double p_init[SIZE];
  double p_tran[SIZE][SIZE]; //p_tran[a][b] = p(b|a)
    
  double *qc_init, (*qc_trans)[SIZE];
};

#endif /* QUALDIST_H_ */
