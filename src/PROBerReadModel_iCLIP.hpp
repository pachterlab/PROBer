#ifndef PROBERREADMODEL_ICLIP
#define PROBERREADMODEL_ICLIP

#include<cassert>
#include<string>

#include "utils.h"
#include "RefSeq.hpp"

#include "SEQstring.hpp"
#include "QUALstring.hpp"
#include "CIGARstring.hpp"
#include "MDstring.hpp"

#include "BamAlignment.hpp"
#include "AlignmentGroup.hpp"

#include "SequencingModel.hpp"

class PROBerReadModel_iCLIP {
public:
  // default constructor
  PROBerReadModel_iCLIP();
  
  /*
    @function   constructor
    @param   model_type   0, SE, no qual; 1, SE qual; 2, PE, no qual; 3 PE, qual
    @param   max_len      maximum read length
   */
  PROBerReadModel_iCLIP(int model_type, int max_len = -1);

  /*
    @function   destructor
  */
  ~PROBerReadModel_iCLIP();

  /*
    @function   get model_type
  */
  int getModelType() const { return model_type; }

  void update(AlignmentGroup& ag);
  void finish() { seqmodel->finish(); }
  void calcProbs(AlignmentGroup& ag, double* conprbs);
  
  void read(const char* modelF);
  void write(const char* modelF);
private:
  int model_type;
  SequencingModel* seqmodel;

  // we can only put these variables here provided this piece of code will not run in parallel
  SEQstring seq;
  QUALstring qual;
  CIGARstring cigar;
  MDstring mdstr;

  RefSeq refseq;
};

inline void PROBerReadModel_iCLIP::update(AlignmentGroup& ag) {
  int size = ag.size();
  BamAlignment *ba = NULL;
  char dir;

  assert(ag.getSEQ(seq));
  if (model_type & 1) assert(ag.getQUAL(qual));
  for (int i = 0; i < size; ++i) {
    ba = ag.getAlignment(i);
    dir = ba->getMateDir();
    assert(ba->getCIGAR(cigar));
    assert(ba->getMD(mdstr));
    refseq.setUp(dir, cigar, mdstr, seq);
    seqmodel->update(1.0, dir, 0, &refseq, &cigar, &seq, ((model_type & 1) ? &qual : NULL));
  }

  if (model_type >= 2) {
    assert(ag.getSEQ(seq, 2));
    if (model_type & 1) assert(ag.getQUAL(qual, 2));
    for (int i = 0; i < size; ++i) {
      ba = ag.getAlignment(i);
      dir = ba->getMateDir(2);
      assert(ba->getCIGAR(cigar, 2));
      assert(ba->getMD(mdstr, 2));
      refseq.setUp(dir, cigar, mdstr, seq);
      seqmodel->update(1.0, dir, 0, &refseq, &cigar, &seq, ((model_type & 1) ? &qual : NULL));
    }
  }
}

inline void PROBerReadModel_iCLIP::calcProbs(AlignmentGroup& ag, double* conprbs) {
  int size = ag.size();
  BamAlignment *ba = NULL;
  char dir;

  assert(ag.getSEQ(seq));
  if (model_type & 1) assert(ag.getQUAL(qual));
  for (int i = 0; i < size; ++i) {
    ba = ag.getAlignment(i);
    dir = ba->getMateDir();
    assert(ba->getCIGAR(cigar));
    assert(ba->getMD(mdstr));
    refseq.setUp(dir, cigar, mdstr, seq);
    conprbs[i] = seqmodel->getProb(dir, 0, &refseq, &cigar, &seq, ((model_type & 1) ? &qual : NULL));
  }

  if (model_type >= 2) {
    assert(ag.getSEQ(seq, 2));
    if (model_type & 1) assert(ag.getQUAL(qual, 2));
    for (int i = 0; i < size; ++i) {
      ba = ag.getAlignment(i);
      dir = ba->getMateDir(2);
      assert(ba->getCIGAR(cigar, 2));
      assert(ba->getMD(mdstr, 2));
      refseq.setUp(dir, cigar, mdstr, seq);
      conprbs[i] *= seqmodel->getProb(dir, 0, &refseq, &cigar, &seq, ((model_type & 1) ? &qual : NULL));
    }
  }

  double sum = 0.0;
  for (int i = 0; i < size; ++i) sum += conprbs[i];

  //assert(sum > 0.0);
  if (sum <= 0.0) sum = 1.0;
  
  for (int i = 0; i < size; ++i) conprbs[i] /= sum;  
}

#endif
