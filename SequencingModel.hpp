#ifndef SEQUENCINGMODEL_H_
#define SEQUENCINGMODEL_H_

#include<cmath>
#include<cassert>
#include<fstream>

#include "RefSeq.hpp"
#include "CIGARstring.hpp"
#include "SEQstring.hpp"
#include "QUALstring.hpp"

#include "Markov.hpp"
#include "Profile.hpp"
#include "QProfile.hpp"

class SequencingModel {
public:
  SequencingModel(bool = true, int = 1000);
  ~SequencingModel();

  double getProb(int, const RefSeq&, const CIGARstring*, const SEQstring*, const QUALstring* = NULL);
  void update(double, int, const RefSeq& refseq, const CIGARstring* cigar, const SEQstring*, const QUALstring* = NULL);

  void init();
  void collect(const SequencingModel*);
  void finish();

  void read(std::ifstream&);
  void write(std::ofstream&);

  void simulate(Sampler*, int, int, const RefSeq&, const std::string&, std::string&, std::string&);

  void startSimulation();
  void finishSimulation();

private:
  bool hasQual; // if quality scores are available

  Markov *markov;
  Profile *profile;
  QProfile *qprofile;

  void push(std::string& cigar, char opchr, int oplen) {
    while (oplen > 0) { cigar.push_back(oplen % 10 + '0'); oplen /= 10; }
    cigar.push_back(opchr);
  }
};

// We assume the direction is set up 
inline double SequencingModel::getProb(int pos, const RefSeq& refseq, const CIGARstring* cigar, const SEQstring* seq, const QUALstring* qual) {
  double prob = 1.0;
  int len = cigar->getLen();
  int readpos = 0;

  char opchr, last_opchr = 0;
  int oplen;

  for (int i = 0; i < len; ++i) {
    opchr = cigar->opchrAt(i);
    oplen = cigar->oplenAt(i);

    // Markov model probabilities
    prob *= (last_opchr == 0 ? markov->getProb(opchr) : markov->getProb(last_opchr, opchr));
    if (oplen > 1) prob *= pow(markov->getProb(opchr, opchr), oplen - 1);

    if (opchr == 'M' || opchr == '=' || opchr == 'X') {
      for (int j = 0; j < oplen; ++j) {
	int ref_base = refseq.baseCodeAt(pos);
	int read_base = seq->baseCodeAt(readpos);
	prob *= (hasQual ? qprofile->getProb(qual->qualAt(readpos), ref_base, read_base) : profile->getProb(readpos, ref_base, read_base));
	++pos; ++readpos;
      }
    }
    else if (opchr == 'I') {
      for (int j = 0; j < oplen; ++j) {
	prob *= markov->getIBaseProb(seq->baseCodeAt(readpos));
	++readpos;
      }
    }
    else {
      assert(opchr == 'D');
      pos += oplen;
    }

    last_opchr = opchr;
  }

  return prob;
}

// We assume the direction is set up 
inline void SequencingModel::update(double frac, int pos, const RefSeq& refseq, const CIGARstring* cigar, const SEQstring* seq, const QUALstring* qual) {
  int len = cigar->getLen();
  int readpos = 0;

  char opchr, last_opchr = 0;
  int oplen;

  for (int i = 0; i < len; ++i) {
    opchr = cigar->opchrAt(i);
    oplen = cigar->oplenAt(i);

    // Markov model probabilities
    if (last_opchr == 0) markov->update(opchr, frac);
    else markov->update(last_opchr, opchr, frac);

    if (oplen > 1) markov->update(opchr, opchr, frac * (oplen - 1));

    if (opchr == 'M' || opchr == '=' || opchr == 'X') {
      for (int j = 0; j < oplen; ++j) {
	int ref_base = refseq.baseCodeAt(pos);
	int read_base = seq->baseCodeAt(readpos);
	if (hasQual) qprofile->update(qual->qualAt(readpos), ref_base, read_base, frac);
	else profile->update(readpos, ref_base, read_base, frac);
	++pos; ++readpos;
      }
    }
    else if (opchr == 'I') {
      for (int j = 0; j < oplen; ++j) {
	markov->updateIBase(seq->baseCodeAt(readpos), frac);
	++readpos;
      }
    }
    else {
      assert(opchr == 'D');
      pos += oplen;
    }

    last_opchr = opchr;
  }
}

// len is the length of the read. 
// qual is a string of (SANGER score + 33)
// Here we assume that the insert_length is approximiated as fragment length based on the assumpiton that insertion/deletion are not common
// If the end of sequence is reached, only insertions are generated
// For refseq, we assume its direction is already set up
inline void SequencingModel::simulate(Sampler *sampler, int len, int pos, const RefSeq& refseq, const std::string& qual, std::string& cigar, std::string& seq) {
  int readpos = 0;
  char opchr, last_opchr;
  int oplen;
  int totLen = refseq.getTotLen();

  cigar.clear();
  seq.clear();

  last_opchr = 0; oplen = 0;
  while (readpos < len) {
    if (pos >= totLen) opchr = 'I';
    else if (last_opchr == 0) opchr = markov->simulate(sampler);
    else opchr = markov->simulate(sampler, last_opchr);
    
    if (opchr == 'M') {
      seq.push_back(hasQual ? qprofile->simulate(sampler, qual[readpos] - 33, refseq.baseCodeAt(pos)) : profile->simulate(sampler, readpos, refseq.baseCodeAt(pos)));
      ++readpos; ++pos;
    }
    else if (opchr == 'I') {
      seq.push_back(markov->simulateIBase(sampler));
      ++readpos;
    }
    else {
      assert(opchr == 'D');
      ++pos;
    }

    if (last_opchr != opchr) { 
      if (last_opchr != 0) push(cigar, last_opchr, oplen);
      last_opchr = opchr;
      oplen = 0;
    }
    ++oplen;
  }
  if (last_opchr != 0) push(cigar, last_opchr, oplen);
}

#endif
