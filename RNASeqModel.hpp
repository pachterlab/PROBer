#ifndef RNASEQMODEL_H_
#define RNASEQMODEL_H_

#include<cassert>
#include<string>
#include<fstream>
#include<algorithm>

#include "utils.h"
#include "sampling.hpp"

#include "RefSeq.hpp"
#include "Refs.hpp"
#include "Transcripts.hpp"

#include "SEQstring.hpp"
#include "QUALstring.hpp"
#include "CIGARstring.hpp"

#include "BamAlignment.hpp"
#include "AlignmentGroup.hpp"

#include "Orientation.hpp"
#include "LenDist.hpp"
#include "RSPD.hpp"
#include "SequencingModel.hpp"
#include "NoiseProfile.hpp"

class RNASeqModel {
public:
  /*
    @param     refs   a pointer to the reference sequences
    @comment:  Can be used to load a learned model file
   */
  RNASeqModel(Refs* refs = NULL);
  
  /*
    @param  paramsF   File contains initial parameters for each component of the model
       paramsF has the following format:

       model_type
       minL maxL mean sd  // fragment length distribution
       mate1_minL mate1_maxL mate2_minL mate2_maxL
       probF // for orientation
       imdName              

    @param  option    0, preprocessing; 1, master; 2, slave
    @param  refs   a pointer to reference sequences
   */
  RNASeqModel(const char* paramsF, int option, Refs* refs = NULL);

  ~RNASeqModel();

  void setE2IArray(const int* e2i) {
    this->e2i = e2i;
  }

  void setSampler(Sampler *sampler) { 
    this->sampler = sampler;
  }

  void update_preprocess(AlignmentGroup& ag, bool isAligned);
  
  /*
    @output outF     The first section is MLD1. The second section is MLD2 if reads are paired-end. 
                     The next section is QualDist if reads contain quality scores. 
		     The last section is for NoiseProfile
   */
  void write_preprocess(const char* imdName);
  void read_preprocess(const char* imdName);

  void init();
  void collect(RNASeqModel* o);
  void finish();

  void read(const char* modelF);
  void write(const char* modelF);

  void simulate(READ_INT_TYPE rid, int& sid, std::ofstream* out1, std::ofstream* out2 = NULL);

  void startSimulation(Sampler* sampler, const std::vector<double>& theta);
  void finishSimulation();

  int getModelType() const { return model_type; }

  // Polish theta for polyA tail effects
  void polishTheta(std::vector<double>& theta);

  // Calculate expected effective lengths
  // This normalization is not "consistent" with the generative model based on theta, but maybe a good approximation
  void calcEEL(std::vector<double>& eel);

  // Calculate a better normalizer, it should be the expected number of fragments generated given a length
  // However, here we just use the length / fragment length as an approximator
  void calcNormalizer(std::vector<double>& normalizers);

private:
  static const double MINEEL = 1.0; // Minimum eel allowed

  int model_type; // 0, SE, no Qual; 1, SE, Qual; 2, PE, no Qual; 3, PE, Qual

  SEQstring seq;
  QUALstring qual;
  CIGARstring cigar;

  Refs *refs;
  Transcripts *transcripts;

  Orientation *ori;
  LenDist *fld, *mld1, *mld2; // fld, fragment length distribution; mld1, mate length distribution 1; mld2, mate length distribution 2.
  RSPD *rspd;
  QualDist *qd;
  SequencingModel *seqmodel;
  NoiseProfile *npro;
  
  Sampler *sampler;
  double *theta_cdf; // for simulation

  int M;
  double *mw; // for masking, only used when poly(A) tails are added

  const int *e2i; // external id to internal id array 

  int range_approx;

  void calcMW();
};

inline void RNASeqModel::update_preprocess(AlignmentGroup& ag, bool isAligned) {
  // Update MLDs
  int len = ag.getSeqLength(1);
  mld1->update(len, 1.0);
  if (!isAligned) mld1->updateC(len);
  if (model_type >= 2) {
    len = ag.getSeqLength(2);
    mld2->update(len, 1.0);
    if (!isAligned) mld2->updateC(len);
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

inline void RNASeqModel::simulate(READ_INT_TYPE rid, int& sid, std::ofstream* out1, std::ofstream* out2) {
  char dir;
  int pos, m2pos, effL;
  int fragLen, mateL1, mateL2;
  std::string qual1, qual2, cigar1, cigar2, readseq1, readseq2;
  
  // Simulate reads
  sid = sampler->sample(theta_cdf, M + 1);
  if (sid == 0) {
    dir = pos = fragLen = 0;
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
    RefSeq &ref = refs->getRef(sid);
    dir = ori->simulate(sampler);
    ref.setDir(dir);

    fragLen = (fld != NULL ? fld->simulateAdjusted(sampler, ref.getTotLen()) : mld1->simulateAdjusted(sampler, ref.getTotLen()));
    assert(fragLen > 0);
    effL = std::min(ref.getFullLen(), ref.getTotLen() - fragLen + 1);
    pos = rspd->simulate(sampler, effL);
    assert(pos >= 0);
    if (dir == '-') pos = ref.getTotLen() - pos - fragLen;
    
    if (fld == NULL) { mateL1 = fragLen; fragLen = 0; }
    else { mateL1 = mld1->simulateAdjusted(sampler, fragLen); assert(mateL1 > 0); }
    if (model_type & 1) qd->simulate(sampler, mateL1, qual1);
    seqmodel->simulate(sampler, mateL1, pos, ref, qual1, cigar1, readseq1);

    if (model_type >= 2) {
      mateL2 = mld2->simulateAdjusted(sampler, fragLen); 
      assert(mateL2 > 0);
      if (model_type & 1) qd->simulate(sampler, mateL2, qual2);
      ref.setDir(dir == '+' ? '-' : '+');
      m2pos = ref.getTotLen() - pos - fragLen;
      seqmodel->simulate(sampler, mateL2, m2pos, ref, qual2, cigar2, readseq2);
    }
  }

  // Output reads
  (*out1)<< ((model_type & 1) ? '@' : '>') << rid<< '_'<< sid<< '_'<< dir<< '_'<< pos<< '_'<< fragLen<< '_'<< cigar1<< (model_type >= 2 ? "/1" : "") << std::endl;
  (*out1)<< readseq1<< std::endl;
  if (model_type & 1) (*out1)<< '+'<< std::endl<< qual1<< std::endl;

  if (model_type >= 2) {
    (*out2)<< ((model_type & 1) ? '@' : '>') << rid<< '_'<< sid<< '_'<< dir<< '_'<< pos<< '_'<< fragLen<< '_'<< cigar2<< "/2"<< std::endl;
    (*out2)<< readseq2<< std::endl;
    if (model_type & 1) (*out2)<< '+'<< std::endl<< qual2<< std::endl;
  }
}

#endif
