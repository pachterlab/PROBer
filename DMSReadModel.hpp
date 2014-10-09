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

  
  void init();
  void collect(DMSReadModel* o);
  void finish();

  void read(const char* modelF);
  void write(const char* modelF);

  void simulate(READ_INT_TYPE rid, int tid, int pos, std::ofstream* out1, int fragment_length = 0, std::ofstream* out2 = NULL);

  void startSimulation();
  void finishSimulation();

private:
  int model_type; // 0, SE, no Qual; 1, SE, Qual; 2, PE, no Qual; 3, PE, Qual

  SEQstring seq;
  QUALstring qual;
  CIGARstring cigar;

  MateLenDist *mld1, *mld2; // mld1, mate length distribution 1; mld2, mate length distribution 2.
  QualDist *qd;
  SequencingModel *seqmodel;
  NoiseProfile *npro;

  int max_len; // maximum mate length

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
    ag.getQUAL(qual); qd->update(qual);
    if (model_type == 3) {
      ag.getQUAL(qual, 2); qd->update(qual);
    }
  }
  
  // Update NoiseProfile
  if (!isAligned) {
    ag.getSEQ(seq); npro->updateC(seq);
    if (model_type >= 2) {
      ag.getSEQ(seq, 2); npro->updateC(seq);
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
    seqmodel->simulate(sampler, mateL1, pos, ref, qual1, cigar1, readseq1);

    if (model_type >= 2) {
      mateL2 = mld2->simulate(sampler, fragment_length); 
      if (model_type & 1) qd->simulate(sampler, mateL2, qual2);
      ref.setDir('-');
      m2pos = ref.getTotLen() - pos - fragment_length;
      seqmodel->simulate(sampler, mateL2, m2pos, ref, qual2, cigar2, readseq2);
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
