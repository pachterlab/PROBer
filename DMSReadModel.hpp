#ifndef DMSREADMODEL_H_
#define DMSREADMODEL_H_

#include<cassert>
#include<string>
#include<fstream>
#include<algorithm>

#include "utils.h"
#include "sampling.hpp"

#include "RefSeq.hpp"
#include "Refs.hpp"

#include "SEQstring.hpp"
#include "QUALstring.hpp"
#include "CIGARstring.hpp"

#include "BamAlignment.hpp"
#include "AlignmentGroup.hpp"

#include "MateLenDist.hpp"
#include "SequencingModel.hpp"
#include "NoiseProfile.hpp"
#include "QualDist.hpp"

#include "InMemoryStructs.hpp"

class DMSReadModel {
public:

  /*
    @param   model_type   0, SE, no qual; 1, SE qual; 2, PE, no qual; 3 PE, qual
    @param   refs         a pointer to reference sequences
    @comment: Master thread for learning model parameters
   */
  DMSReadModel(int model_type, Refs* refs);

  /*
    @param   master_model   The master model for learning parameters
    @comment: Constructor function for slave threads
   */
  DMSReadModel(DMSReadModel* master_model);

  /*
    @param     refs     a pointer to the reference sequences
    @param     sampler  a pointer to a sampler for simulation
    @comment:  Used for simulation only, read learned parameters from files 
   */
  DMSReadModel(Refs* refs, Sampler* sampler);
  
  ~DMSReadModel();

  /*
    @func   get model type
   */
  int getModelType() const { return model_type; }

  void update_preprocess(AlignmentGroup& ag, bool isAligned);

  void finish_preprocess();

  /*
    @param   ag_in_mem   an in-memory alignment group, recorded information necessary for EM iteration
    @param   aligns      a pointer to all alignments in ag_in_mem
    @param   ag   an alignment group, which contains the read sequence etc.
    @func   set conditional probabilities to an alignments of a read
   */
  void setConProbs(InMemAlignG* ag_in_mem, InMemAlign* aligns, AlignmentGroup& ag);

  /*
    @param   ag_in_mem
    @param   aligns
    @param   ag
    @param   noise_frac   fractional weight at noise transcript
   */
  void update(InMemAlignG* ag_in_mem, InMemAlign* aligns, AlignmentGroup& ag, double noise_frac);

  /*
    @return  partial log-likelihood for unalignable reads
   */
  double calcLogP() {
    return loglik + npro->calcLogP();
  }

  void init();
  void collect(DMSReadModel* o);
  void finish();

  void read(const char* modelF);
  void write(const char* modelF);

  void simulate(READ_INT_TYPE rid, int tid, int pos, int fragment_length, std::ofstream* out1, std::ofstream* out2 = NULL);

  void startSimulation();
  void finishSimulation();

private:
  int model_type; // 0, SE, no Qual; 1, SE, Qual; 2, PE, no Qual; 3, PE, Qual

  MateLenDist *mld1, *mld2; // mld1, mate length distribution 1; mld2, mate length distribution 2.
  QualDist *qd;
  SequencingModel *seqmodel;
  NoiseProfile *npro;

  int max_len; // maximum mate length
  double loglik; // partial log-likelihood for unaligned reads

  Refs *refs;  
  Sampler *sampler;
};

inline void DMSReadModel::update_preprocess(AlignmentGroup& ag, bool isAligned) {
  // Update MLDs
  int len = ag.getSeqLength(1);
  mld1->update(len, !isAligned);
  if (model_type >= 2) {
    len = ag.getSeqLength(2);
    mld2->update(len, !isAligned);
  }
  
  // Updae QualDist
  if (model_type & 1) {
    QUALstring qual;
    ag.getQUAL(qual); qd->update(qual);
    if (model_type == 3) {
      ag.getQUAL(qual, 2); qd->update(qual);
    }
  }
  
  // Update NoiseProfile
  if (!isAligned) {
    SEQstring seq;
    ag.getSEQ(seq); npro->updateC(seq);
    if (model_type >= 2) {
      ag.getSEQ(seq, 2); npro->updateC(seq);
    }
  }
}

inline void DMSReadModel::setConProbs(InMemAlignG* ag_in_mem, InMemAlign* aligns, AlignmentGroup& ag) {
  int seqlen;
  SEQstring seq;
  QUALstring qual;
  CIGARstring cigar;

  // Get read sequences and quality scores
  assert(ag.getSEQ(seq));
  seqlen = ag.getSeqLength();
  if (model_type & 1) assert(ag.getQUAL(qual));
  // set noise probability    
  ag_in_mem->noise_conprb = mld1->getProb(seqlen) * npro->getProb(seq);
  // set alignment probabilities
  for (int i = 0; i < ag_in_mem->size; ++i) if (aligns[i].conprb != -1.0) {
    RefSeq &refseq = refs->getRef(aligns[i].tid);
    assert(ag.getAlignment(i)->getCIGAR(cigar));
    aligns[i].conprb = (aligns[i].fragment_length > 0 ? mld1->getProb(seqlen, aligns[i].fragment_length) : mld1->getProb(seqlen)) * \
      seqmodel->getProb('+', aligns[i].pos, refseq, &cigar, &seq, ((model_type & 1) ? &qual : NULL));
  }
  
  if (model_type >= 2) {
    // paired-end reads
    assert(ag.getSEQ(seq, 2));
    seqlen = ag.getSeqLength(2);
    if (model_type & 1) assert(ag.getQUAL(qual, 2));
    ag_in_mem->noise_conprb *= mld2->getProb(seqlen) * npro->getProb(seq);
    for (int i = 0; i < ag_in_mem->size; ++i) if (aligns[i].conprb != -1.0) {
      RefSeq &refseq = refs->getRef(aligns[i].tid);
      assert(ag.getAlignment(i)->getCIGAR(cigar, 2));
      assert(aligns[i].fragment_length > 0);
      aligns[i].conprb *= mld2->getProb(seqlen, aligns[i].fragment_length) * \
	seqmodel->getProb('-', refseq.getTotLen() - aligns[i].pos - aligns[i].fragment_length, refseq, &cigar, &seq, ((model_type & 1) ? &qual : NULL));
    }
  }
}

inline void DMSReadModel::update(InMemAlignG* ag_in_mem, InMemAlign* aligns, AlignmentGroup& ag, double noise_frac) {
  SEQstring seq;
  QUALstring qual;
  CIGARstring cigar;

  assert(ag.getSEQ(seq));
  if (model_type & 1) assert(ag.getQUAL(qual));
  // update noise prob
  npro->update(seq, noise_frac);
  // update alignment probs
  for (int i = 0; i < ag_in_mem->size; ++i) if (aligns[i].frac > 0.0) {
    RefSeq &refseq = refs->getRef(aligns[i].tid);
    assert(ag.getAlignment(i)->getCIGAR(cigar));
    seqmodel->update(aligns[i].frac, '+', aligns[i].pos, refseq, &cigar, &seq, ((model_type & 1) ? &qual : NULL));
  }

  if (model_type >= 2) {
    // paired-end reads
    assert(ag.getSEQ(seq, 2));
    if (model_type & 1) assert(ag.getQUAL(qual, 2));
    // update noise prob
    npro->update(seq, noise_frac);
    // update alignment probs
    for (int i = 0; i < ag_in_mem->size; ++i) if (aligns[i].frac > 0.0) {
      RefSeq &refseq = refs->getRef(aligns[i].tid);
      assert(ag.getAlignment(i)->getCIGAR(cigar, 2));
      assert(aligns[i].fragment_length > 0);
      seqmodel->update(aligns[i].frac, '-', refseq.getTotLen() - aligns[i].pos - aligns[i].fragment_length, refseq, &cigar, &seq, ((model_type & 1) ? &qual : NULL));
    }
  }
} 

inline void DMSReadModel::simulate(READ_INT_TYPE rid, int tid, int pos, int fragment_length, std::ofstream* out1, std::ofstream* out2) {
  int m2pos;
  int mateL1, mateL2;
  std::string qual1, qual2, cigar1, cigar2, readseq1, readseq2;
  
  // Simulate reads
  if (tid == 0) {
    cigar1 = cigar2 = "*";

    mateL1 = mld1->simulate(sampler);
    if (model_type & 1) qd->simulate(sampler, mateL1, qual1);
    npro->simulate(sampler, mateL1, readseq1);
    
    if (model_type >= 2) {
      mateL2 = mld2->simulate(sampler);
      if (model_type & 1) qd->simulate(sampler, mateL2, qual2);
      npro->simulate(sampler, mateL2, readseq2);
    }
  }
  else {
    RefSeq &ref = refs->getRef(tid);
    mateL1 = mld1->simulate(sampler, fragment_length);
    if (model_type & 1) qd->simulate(sampler, mateL1, qual1);
    seqmodel->simulate(sampler, mateL1, '+', pos, ref, qual1, cigar1, readseq1);

    if (model_type >= 2) {
      mateL2 = mld2->simulate(sampler, fragment_length); 
      if (model_type & 1) qd->simulate(sampler, mateL2, qual2);
      m2pos = ref.getTotLen() - pos - fragment_length;
      seqmodel->simulate(sampler, mateL2, '-', m2pos, ref, qual2, cigar2, readseq2);
    }
  }

  // Output reads
  (*out1)<< ((model_type & 1) ? '@' : '>') << rid<< '_'<< tid<< '_'<< pos<< '_'<< fragment_length<< '_'<< cigar1<< (model_type >= 2 ? "/1" : "") << std::endl;
  (*out1)<< readseq1<< std::endl;
  if (model_type & 1) (*out1)<< '+'<< std::endl<< qual1<< std::endl;

  if (model_type >= 2) {
    (*out2)<< ((model_type & 1) ? '@' : '>') << rid<< '_'<< tid<< '_'<< pos<< '_'<< fragment_length<< '_'<< cigar2<< "/2"<< std::endl;
    (*out2)<< readseq2<< std::endl;
    if (model_type & 1) (*out2)<< '+'<< std::endl<< qual2<< std::endl;
  }
}

#endif
