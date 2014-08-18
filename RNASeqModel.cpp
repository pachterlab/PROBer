#include<cstdio>
#include<cstring>
#include<string>
#include<fstream>
#include<vector>
#include<algorithm>

#include "utils.h"
#include "my_assert.h"

#include "RefSeq.hpp"
#include "Refs.hpp"
#include "Transcripts.hpp"

#include "BamAlignment.hpp"
#include "AlignmentGroup.hpp"

#include "Orientation.hpp"
#include "LenDist.hpp"
#include "RSPD.hpp"
#include "SequencingModel.hpp"
#include "NoiseProfile.hpp"

#include "RNASeqModel.hpp"

RNASeqModel::RNASeqModel(Refs* refs) {
  model_type = -1; // not initilaized yet

  ori = NULL;
  fld = mld1 = mld2 = NULL;
  rspd = NULL;
  qd = NULL;
  seqmodel = NULL;
  npro = NULL;

  sampler = NULL;
  theta_cdf = NULL;

  M = 0; mw = NULL;
  this->refs = refs;
  if (refs != NULL) M = refs->getM();

  e2i = NULL;

  range_approx = 201;
}

RNASeqModel::RNASeqModel(const char* paramsF, int option, Refs* refs) {
  int minL, maxL;
  std::ifstream fin(paramsF);

  M = 0; mw = NULL;
  this->refs = refs;
  if (refs != NULL) M = refs->getM();

  assert(fin.is_open());  
  assert(fin>> model_type);
  assert(fin>> minL>> maxL);

  if (option == 0) {
    // pre-estimate a set of parameters when parsing the alignment file
    mld1 = new LenDist(minL, maxL, true);
    if (type >= 2) mld2 = new LenDist(minL, maxL, true);
    if (model_type == 1 || model_type == 3) qd = new QualDist();
    npro = new NoiseProfile(true);
  }
  else {
    double mean, sd, probF;
    int mate1_minL, mate1_maxL, mate2_minL, mate2_maxL;
    char imdName[STRLEN];

    assert(fin>> mean>> sd);
    // If SE reads, mate2_minL = mate2_maxL = 0
    assert(fin>> mate1_minL>> mate1_maxL>> mate2_minL>> mate2_maxL);
    assert(fin>> probF);
    assert(fin>> imdName);
    assert(fin>> range_approx);

    // Components shared by master and slave threads
    ori = new Orientation(probF);
    if (model_type >= 2) fld = new LenDist(minL, maxL);
    seqmodel = new SequencingModel((model_type == 1 || model_type == 3), std::max(mate1_maxL, mate2_maxL));
    npro = new NoiseProfile(option == 1 ? true : false);

    // Components specific to master thread
    if (option == 1) {
      if (model_type < 2 && mean > 0.0) {
	fld = new LenDist(minL, maxL);
	fld->setAsNormal(mean, sd);
      }
      mld1 = new LenDist(mate1_minL, mate1_maxL);
      if (model_type >= 2) mld2 = new LenDist(mate2_minL, mate2_maxL);
      rspd = new RSPD();
      qd = new QualDist();

      // Read pre-estimated parameters
      read_preprocess(imdName);

      if (refs->hasPolyA()) {
	mw = new double[M + 1];
	calcMW();
      }
    }
  }

  fin.close();

  e2i = NULL;
}

RNASeqModel::~RNASeqModel() {
  refs = NULL;
  transcripts = NULL;

  if (ori != NULL) delete ori;
  if (fld != NULL) delete fld;
  if (mld1 != NULL) delete mld1;
  if (mld2 != NULL) delete mld2;
  if (rspd != NULL) delete rspd;
  if (qd != NULL) delete qd;
  if (seqmodel != NULL) delete seqmodel;
  if (npro != NULL) delete npro;

  sampler = NULL;

  if (mw != NULL) delete[] mw;
}

void RNASeqModel::write_preprocess(const char* imdName) {
  char tmpModelF[STRLEN], tmpLogProbF[STRLEN];
  sprintf(tmpModelF, "%s.tmp_model", imdName);
  sprintf(tmpLogProbF, "%s.tmp_logprob", imdName);
  
  std::ofstream fout(tmpModelF);
  double logp_mld1 = 0.0, logp_mld2 = 0.0, logp_qd = 0.0;
  
  assert(fout.is_open());

  // Write out mate length distributions
  mld1->finish();
  mld1->write(fout);
  logp_mld1 = mld1->calcLogP();

  if (model_type >= 2) {
    mld2->finish();
    mld2->write(fout);
    logp_mld2 = mld2->calcLogP();
  }

  // Write out QualDist
  if (model_type & 1) {
    logp_qd = qd->finish();
    qd->write(fout);
  }

  // Write out NoiseProfile
  npro->writeC(fout);
  fout.close();

  // Write down the pre-caluclated log probabilities
  fout.open(tmpLogProbF);
  fout.precision(10);
  fout.setf(0, std::ios::floatfield);
  fout<< logp_mld1<< logp_mld2<< logp_qd<< std::endl;
  fout.close();
}

void RNASeqModel::read_preprocess(const char* imdName) {
  char tmpModelF[STRLEN];
  sprintf(tmpModelF, "%s.tmp_model", imdName);

  std::ifstream fin(tmpModelF);
  assert(fin.is_open());
  
  // Read mate length distributions
  mld1->read(fin);
  if (model_type >= 2) mld2->read(fin);
  
  // Read QualDist
  if (model_type & 1) qd->read(fin);
  
  // Read NoiseProfile
  npro->readC(fin);

  fin.close();
}

void RNASeqModel::init() {
  ori->init();
  if (model_type >= 2) fld->init();
  if (rspd->estimateRSPD()) rspd->init();
  seqmodel->init();
  npro->init();
}

void RNASeqModel::collect(RNASeqModel* o) {
  ori->collect(o->ori);
  if (model_type >= 2) fld->collect(o->fld);
  if (rspd->estimateRSPD()) rspd->collect(o->rspd);
  seqmodel->collect(o->seqmodel);
  npro->collect(o->npro);
}

void RNASeqModel::finish() {
  ori->finish();
  if (model_type >= 2) fld->finish();
  if (rpsd->estimateRSPD()) rspd->finish();
  seqmodel->finish();
  npro->finish();
  
  if ((mw != NULL) && (model_type >= 2 || rspd->estimateRSPD())) calcMW();
}

void RNASeqModel::read(const char* modelF) {
  std::string line;
  std::ifstream fin(modelF);
  assert(fin.is_open());

  // Read model type
  while (getline(fin, line)) {
    if (line.substr(0, 11) == "#Model Type") break;
  }
  assert(fin.good());
  assert(fin>> model_type);
  getline(fin, line);

  // Read orientation
  ori = new Orientation();
  ori->read(fin);

  // Read fragment length distribution
  while (getline(fin, line)) {
    if (line == "#Fragment length distribution") break;
  }
  assert(fin.good());
  assert(getline(fin, line));
  if (line == "#yes") {
    fld = new LenDist();
    fld->read(fin);
  }
  else 
    assert((line == "#no") && (model_type < 2));
    
  // Read mate length distribution 1
  mld1 = new LenDist();
  mld1->read(fin);

  // Read mate length distribution 2
  if (model_type >= 2) {
    mld2 = new LenDist();
    mld2->read(fin);
  }

  // Read RSPD
  rspd = new RSPD();
  rspd->read(fin);

  // Read QualDist
  if (model_type & 1) {
    qd = new QualDist();
    qd->read(fin);
  }

  // Read sequencing model
  seqmodel = new SequencingModel((model_type & 1));
  seqmodel->read(fin);

  // Read noise profile model
  npro = new NoiseProfile();
  npro->read(fin);

  // Read mw values, this one is optional
  int value;
  if (fin>> value) {
    if (M == 0) M = value;
    if (M == value) {
      mw = new double[M + 1];
      for (int i = 0; i <= M; ++i) assert(fin>> mw[i]);
    }
  }

  fin.close();
}

void RNASeqModel::write(const char* modelF) {
  std::ofstream fout(modelF);
  assert(fout.is_open());

  // Write model type
  fout<< "#Model Type: 0, SE, no qual; 1, SE, qual; 2, PE, no qual; 3 PE, qual"<< std::endl;
  fout<< model_type<< std::endl<< std::endl;

  // Write orientation
  ori->write(fout);
  
  // Write fragment length distribution
  fout<< "#Fragment length distribution"<< std::endl;
  fout<< (fld == NULL ? "#no" : "#yes")<< std::endl;
  if (fld != NULL) fld->write(fout);

  // Write mate length distribution 1
  fout<< "#Mate length distribution 1"<< std::endl;
  mld1->write(fout);

  // Write mate length distribution 2
  if (model_type >= 2) {
    fout<< "#Mate length distribution 2"<< std::endl;
    mld2->write(fout);
  }

  // Write RSPD
  rspd->write(fout);

  // Write QualDist
  if (model_type & 1) qd->write(fout);

  // Write sequencing model
  seqmodel->write(fout);

  // Write noise profile
  npro->write(fout);

  // Write mw values
  if (mw != NULL) {
    fout<< M<< std::endl;
    fout.precision(10);
    fout.setf(0, std::ios::floatfield);
    for (int i = 0; i < M; ++i) fout<< mw[i]<< '\t';
    fout<< mw[M]<< std::endl;
  }

  fout.close();
}

void RNASeqModel::startSimulation(Sampler* sampler, const std::vector<double>& theta) {
  this->sampler = sampler;
  
  theta_cdf = new double[M + 1];
  theta_cdf[0] = theta[0];
  for (int i = 1; i <= M; ++i) 
    theta_cdf[i] = theta_cdf[i - 1] + theta[i];

  if (model_type & 1) qd->startSimulation();
  seqmodel->startSimulation();
  npro->startSimulation();
}

void RNASeqModel::finishSimulation() {
  delete[] theta_cdf;
  
  if (model_type & 1) qd->finishSimulation();
  seqmodel->finishSimulation();
  npro->finishiSimulation();
}

void RNASeqModel::polishTheta(const std::vector<double>& theta) {
  if (mw == NULL) return;
  double sum = 0.0;
  for (int i = 0; i <= M; ++i) 
    if (isZero(mw[i])) theta[i] = 0.0;
    else {
      theta[i] /= mw[i];
      sum += theta[i];
    }
  assert(!isZero(sum));
  for (int i = 0; i <= M; ++i) theta[i] /= sum;
}

/*
  Expected effective length of a transcript is calculated as:
  EEL := \sigma_{i=minL}^{maxL} pmf[i] * min(fullLen, max(totLen - i + 1, 0))

  For totLen - i + 1 >= fullLen => i <= totLen - fullLen + 1, we have pmf[i] * fullLen term
  Thus for fragment length <= len1, len1 = min(totLen - fullLen + 1, maxL), the sum can be written as cdf[len1] * fullLen
  
  len2 = min(totLen, maxL), which is the maximum fragment this transcript can hold.

  For fragments in the range of [len1 + 1, len2], their sum is \sigma pmf[i] * (totLen - i + 1) = (cdf[len2] - cdf[len1]) * (totLen + 1) - (\simga pmf[i] * i).
  clen[i] = \sigma_{j=1}^{i} pmf[j+minL - 1] * i
 */
void RNASeqModel::calcEEL(std::vector<double>& eel) {
  LenDist *len_dist = NULL;
  double *clen = NULL;
  int minL, maxL;
  int totLen, fullLen, len1, len2;

  len_dist = (fld != NULL ? fld : mld1);
  minL = len_dist->getMinL();
  maxL = len_dist->getMaxL();

  clen = new double[maxL - minL + 2];

  clen[0] = 0.0;
  for (int i = minL; i <= maxL; ++i) {
    int id = i - minL + 1;
    clen[id] = clen[id - 1] + len_dist->getProb(i) * i;
  }
  
  eel.assign(M + 1, 0.0);
  for (int i = 1; i <= M; ++i) {
    totLen = refs->getRef(i).getTotLen();
    fullLen = refs->getRef(i).getFullLen();
    len1 = std::min(totLen - fullLen + 1, maxL);
    len2 = std::min(totLen, maxL);
    
    if (len2 < minL) { eel[i] = 0.0; continue; }
    
    eel[i] = len_dist->getCProb(len1) * fullLen + (len_dist->getSProb(len1 + 1, len2) * (totLen + 1) - (clen[std::max(len2 - minL + 1, 0)] - clen[std::max(len1 - minL + 1, 0)]));
    assert(eel[i] >= 0.0);
    
    if (eel[i] < MINEEL) { eel[i] = 0.0; }
  }

  delete[] clen;
}

// If polyA tail is added, add at most one read
void RNASeqModel::calcNormalizer(std::vector<double>& normalizers) {
  LenDist *len_dist = (fld != NULL ? fld : mld1);
  double meanL = len_dist->calcMean();

  normalizers.assign(M, 0.0);
  for (int i = 1; i <= M; ++i) {
    int fullLen = refs->getRef(i).getFullLen();
    int totLen = refs->getRef(i).getTotLen();
    normalizers[i] = fullLen / meanL;
    if (totLen > fullLen) normalizers[i] = std::min(totLen / meanL, normalizers[i] + 1.0);
  }
}

void RNASeqModel::calcMW() {
  double probF, probR;
  int totLen, fullLen;
  int pos, start, end;
  int fragLen, minL, maxL, upperL;
  double value, value2;
  int start2, minL2, minL_approx, maxL_approx;
  
  probF = ori->getProb('+');
  probR = ori->getProb('-');
  minL = (fld != NULL ? fld->getMinL() : mld1->getMinL());
  maxL = (fld != NULL ? fld->getMaxL() : mld1->getMaxL());
  minL2 = mld1->getMinL();

  mw[0] = 1.0;
  for (int i = 1; i <= M; ++i) {
    mw[i] = 1.0;

    RefSeq& ref = refs->getRef(i);
    totLen = ref.getTotlen();
    fullLen = ref.getFullLen();

    // This transcript is not added a poly(A) tail, skip
    if (fullLen == totLen) continue;
    // This transcript is too short
    if (fullLen <= MASK_LEN) { mw[i] = 0.0; continue; }

    value = 0.0;
    start = fullLen - MASK_LEN;
    end = std::min(fullLen, totLen - minL + 1);    
    for (pos = start; pos < end; ++pos) {
      upperL = std::min(maxL, totLen - pos);
      for (fragLen = minL; fragLen <= upperL; ++fragLen) {
	value += (fld != NULL ? fld->getAdjustedProb(fragLen, totLen) : mld1->getAdjustedProb(fragLen, totLen)) * rspd->getProb(pos, std::min(fullLen, totLen - fragLen + 1));
      }
    }
    
    // SE read, reverse strand
    if (model_type < 2 && fld != NULL) {
      value *= probF;

      // Locate should be done before
      minL_approx = fld->get_len_l();
      maxL_approx = fld->get_len_r();
      if (maxL_approx > totLen) { 
	minL_approx = std::max(minL, minL_approx - (maxL_approx - totLen));
	maxL_approx = totLen;
      }

      value2 = 0.0;
      start2 = start + minL2;
      if (start2 <= totLen) {
	start2 = std::max(0, start2 - maxL_approx);
	end = std::min(fullLen, totLen - minL_approx + 1);
	for (pos = start2; pos < end; ++pos) {
	  bottomL = std::max(minL_approx, start + minL2 - pos);
	  upperL = std::min(maxL_approx, totLen - pos);
	  for (fragLen = bottomL; fragLen <= upperL; ++fragLen) 
	    value2 += fld->getAdjustedProb(fragLen, totLen) * rspd->getProb(pos, std::min(fullLen, totLen - fragLen + 1)) * mld1->getAdjustedCumulativeProb(pos + fragLen - start, fragLen);
	}
      }
      value2 *= probR;
      value += value2;
    }
    
    mw[i] -= value;
    if (isZero(mw[i])) mw[i] = 0.0;
  }
}
