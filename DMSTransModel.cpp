#include<cmath>
#include<cstring>
#include<cassert>
#include<string>
#include<fstream>
#include<algorithm>

#include "utils.h"
#include "sampling.hpp"
#include "DMSTransModel.hpp"

const double DMSTransModel::INF = 1000.0;

int DMSTransModel::primer_length = 6; // default, 6bp
int DMSTransModel::min_frag_len;
int DMSTransModel::max_frag_len;

double DMSTransModel::gamma_init;
double DMSTransModel::beta_init;

double DMSTransModel::base;
double DMSTransModel::cgamma;
double DMSTransModel::dgamma;
double DMSTransModel::cbeta;
double DMSTransModel::dbeta;

int DMSTransModel::min_alloc_len;
bool DMSTransModel::isMAP = true; // default is true

bool DMSTransModel::learning = false; // default is simulation
int DMSTransModel::state = 0;

void DMSTransModel::setGlobalParams(int primer_length, int min_frag_len, int max_frag_len, int init_state) { 
  assert(primer_length <= min_frag_len && min_frag_len <= max_frag_len);
  DMSTransModel::primer_length = primer_length;
  DMSTransModel::min_frag_len = min_frag_len - primer_length;
  DMSTransModel::max_frag_len = max_frag_len - primer_length;
  state = init_state;
}

void DMSTransModel::setLearningRelatedParams(double gamma_init, double beta_init, double base, int read_length, bool isMAP) {
  learning = true;

  DMSTransModel::gamma_init = gamma_init;
  DMSTransModel::beta_init = beta_init;
  DMSTransModel::base = base;
  DMSTransModel::min_alloc_len = (read_length < min_frag_len ? min_frag_len : read_length) - primer_length;
  DMSTransModel::isMAP = isMAP;

  if (isMAP) {
    dgamma = gamma_init * base;
    cgamma = base - dgamma;
    dbeta = beta_init * base;
    cbeta = base - dbeta;
  }
}



DMSTransModel::DMSTransModel(int tid, const std::string& name, int transcript_length) : tid(tid), name(name) {  
  gamma = beta = NULL;
  start = end = NULL;
  dcm = ccm = NULL;
  end_se = NULL;

  logsum = margin_prob = NULL;
  margin_prob2 = NULL;

  start2 = end2 = NULL;
  for (int i = 0; i < 2; ++i) alignmentsArr[i].clear();

  len = efflen = -1; 
  efflen2 = -1;
  N_obs[0] = N_obs[1] = 0.0;
  prob_pass[0] = prob_pass[1] = 1.0; // In case no alignments, unobserved read counts is 0

  delta = 0.0;

  hasSE = false;
  N_se = 0.0;

  cdf_end = NULL;

  if (!learning) return;

  len = transcript_length - primer_length;
  efflen = len - min_frag_len + 1;
  efflen2 = len - min_alloc_len + 1;
  delta = 1.0 / (len + 1.0);
  
  // If a transcript is excluded from analysis, all its gamma/beta values become 0
  gamma = new double[len + 1];
  memset(gamma, 0, sizeof(double) * (len + 1));

  if (getState() > 0) {
    beta = new double[len + 1];
    memset(beta, 0, sizeof(double) * (len + 1));
  }
}

DMSTransModel::~DMSTransModel() {
  if (gamma != NULL) delete[] gamma;
  if (beta != NULL) delete[] beta;

  if (efflen <= 0) return;

  if (start != NULL) delete[] start;
  if (end != NULL) delete[] end;

  // auxiliary arrays
  if (logsum != NULL) delete[] logsum;
  if (margin_prob != NULL) delete[] margin_prob;

  if (margin_prob2 != NULL) delete[] margin_prob2;

  if (dcm != NULL) delete[] dcm;
  if (ccm != NULL) delete[] ccm;

  if (end_se != NULL) delete[] end_se;
}

void DMSTransModel::init() {
  int state = getState();

  // set initial values for EM
  if ((state & 1) == 0) for (int i = 1; i <= len; ++i) gamma[i] = gamma_init;
  if ((state & 1) == 1) for (int i = 1; i <= len; ++i) beta[i] = beta_init;

  if (state < 3) {
    // count vectors for fragment starts and ends
    start = new double[len + 1];
    end = new double[len + 1];
    memset(start, 0, sizeof(double) * (len + 1));
    memset(end, 0, sizeof(double) * (len + 1));
    
    // auxiliary arrays
    logsum = new double[len + 1];
    margin_prob = new double[efflen];
    memset(logsum, 0, sizeof(double) * (len + 1));
    memset(margin_prob, 0, sizeof(double) * efflen);
    
    if (efflen != efflen2 && efflen2 > 0) {
      margin_prob2 = new double[efflen2];
      memset(margin_prob2, 0, sizeof(double) * (efflen2));
    }
    
    if (state == 2) {
      dcm = new double[len + 1];
      ccm = new double[len + 1];
      memset(dcm, 0, sizeof(double) * (len + 1));
      memset(ccm, 0, sizeof(double) * (len + 1));
    }
    
    if (hasSE) {
      end_se = new double[len + 1];
      memset(end_se, 0, sizeof(double) * (len + 1));
    }
  }

  // Calculate auxiliary arrays before the first EM run
  calcAuxiliaryArrays(state & 1);
}

void DMSTransModel::update() {
  int channel = getChannel();
  std::vector<InMemAlign*> &alignments = alignmentsArr[channel];

  HIT_INT_TYPE size = alignments.size();
  double N_se;

  // initialize
  N_obs[channel] = 0.0;
  memset(start, 0, sizeof(double) * (len + 1));
  memset(end, 0, sizeof(double) * (len + 1));

  N_se = 0.0;
  if (hasSE) memset(end_se, 0, sizeof(double) * (len + 1));

  for (HIT_INT_TYPE i = 0; i < size; ++i) {
    end[alignments[i]->pos] += alignments[i]->frac;
    if (alignments[i]->fragment_length > 0) {
      start[alignments[i]->pos + alignments[i]->fragment_length - primer_length] += alignments[i]->frac; 
    }
    else {
      end_se[alignments[i]->pos] += alignments[i]->frac;
      N_se += alignments[i]->frac;
    }
    N_obs[channel] += alignments[i]->frac;
  }

  if (isZero(N_obs[channel])) N_obs[channel] = 0.0; // if N_obs is small, directly set it to 0
}

inline void DMSTransModel::solveQuadratic1(double& beta, double gamma, double dc, double cc) {
  double a = (1.0 - gamma) * (cbeta + cc + dbeta + dc);
  double b = ((cbeta + cc + 2.0 * dbeta + dc) * gamma - (dc + dbeta)) / a;
  double c = (-dbeta * gamma) / a;
  double sqt_delta = sqrt(b * b - 4.0 * c);
  assert(sqt_delta > fabs(b));
  beta = (-b + sqt_delta) / 2.0;
  assert(beta > 0.0 && beta < 1.0);
}

inline void DMSTransModel::solveQuadratic2(double& gamma, double& beta, double dcm, double ccm, double dcp, double ccp) {
  double common_factor = cgamma + ccm - cbeta- dbeta;
  double a = (cbeta + ccp + dbeta + dcp) * common_factor;
  double b = (cbeta + ccp + dbeta) * (dbeta + dgamma + dcm) - common_factor * (dcp + dbeta) + dbeta * dcp;
  double c = - dbeta * (dbeta + dcp + dgamma + dcm);
  double sqt_delta;

  // solve beta
  if (!isZero(fabs(a))) {
    // a != 0
    b /= a; c /= a;
    sqt_delta = sqrt(b * b - 4.0 * c);
    assert(sqt_delta >= 0.0);
    beta = (-b + (a > 0 ? sqt_delta : -sqt_delta)) / 2.0;
  }
  else {
    // a == 0
    beta = - c / b;
  }
  assert(beta > 0.0 && beta < 1.0);

  // calculate gamma given beta
  gamma = (dgamma + dcm) / (cgamma + ccm + dgamma + dcm + (dbeta * (1.0 - beta) / beta - cbeta));
  assert(gamma > 0.0 && gamma < 1.0);
}

void DMSTransModel::EM_step(double N_tot) {
  int channel = getChannel();

  int max_end_i;
  double psum, value;
  
  assert(start2 != NULL && end2 != NULL);

  // What to do if no observed reads
  if (isZero(N_obs[channel])) {
    // force unobserved reads to zero
    switch(getState()) {
    case 0:
      value = (isMAP ? dgamma / (cgamma + dgamma) : 0.0);
      for (int i = 1; i <= len; ++i) gamma[i] = value;
      break;
    case 1:
      value = (isMAP ? dbeta / (cbeta + dbeta) : 0.0);
      for (int i = 1; i <= len; ++i) beta[i] = value;
      break;
    case 2:
      memset(dcm, 0, sizeof(double) * (len + 1));
      memset(ccm, 0, sizeof(double) * (len + 1));
      break;
    case 3:
      if (isMAP) {
	value = dbeta / (cbeta + dbeta);
	for (int i = 1; i <= len; ++i) {
	  gamma[i] = (dgamma + dcm[i]) / (cgamma + ccm[i] + dgamma + dcm[i]); 
	  beta[i] = value;
	}
      }
      else {
	for (int i = 1; i <= len; ++i) {
	  gamma[i] = (dcm[i] > 0.0 ? dcm[i] / (dcm[i] + ccm[i]) : 0.0);
	  beta[i] = 0.0;
	}
      }
      break;
    default: assert(false);
    }
  }
  else {
    //E step, if we have reads that do not know their start positions, infer start from end
    if (!isZero(N_se)) {
      double prev, curr;
      double *mp = NULL;
      
      assert(efflen2 > 0);
      
      mp = (min_alloc_len > min_frag_len ? margin_prob2 : margin_prob);
      curr = (end_se[0] > 0.0 && mp[0] > 0.0) ? end_se[0] / mp[0] : 0.0;
      start[min_alloc_len] += curr;
      for (int i = 1, pos = min_alloc_len + 1; i < efflen2; ++i, ++pos) {
	prev = curr;
	curr = (end_se[i] > 0.0 && mp[i] > 0.0) ? end_se[i] / mp[i] : 0.0;
	max_end_i = (pos - 1) - max_frag_len;
	
	value = prev;
	if (max_end_i >= 0) {
	  value -= ((end_se[max_end_i] > 0.0 && mp[max_end_i] > 0.0) ? end_se[max_end_i] * (exp(logsum[pos - 1] - logsum[max_end_i + min_alloc_len]) / mp[max_end_i]) : 0.0);
	  if (value < 0.0) value = 0.0;
	}
	
	curr += (channel == 0 ? (1.0 - gamma[pos]) : (1.0 - gamma[pos]) * (1.0 - beta[pos])) * value;
	start[pos] += curr;
      }
    }
    
    // E step, calculate hidden reads
    // Calculate end2
    psum = 1.0;
    for (int i = len; i >= 0; --i) {
      if (i < efflen) end2[i] = std::max(psum - exp(logsum[i + min_frag_len] - logsum[i]) * margin_prob[i], 0.0);
      else end2[i] = psum;
      if (i > 0) end2[i] *= (channel == 0 ? gamma[i] : gamma[i] + beta[i] - gamma[i] * beta[i]);
      end2[i] *= delta * N_tot;
      if (i > 0) psum = 1.0 + psum * (channel == 0 ? (1.0 - gamma[i]) : (1.0 - gamma[i]) * (1.0 - beta[i]));
    }
    
    // Calculate start2
    for (int i = 0; i < min_frag_len; ++i) start2[i] = delta * N_tot;
    psum = 1.0;
    for (int i = min_frag_len, pos = 0; i <= len; ++i, ++pos) {
      start2[i] = std::max(1.0 - psum * exp(logsum[i] - logsum[pos]), 0.0);
      start2[i] *= delta * N_tot;
      if (i < len) {
	max_end_i = i - max_frag_len;
	if (max_end_i >= 0) {
	  value = exp(logsum[pos] - logsum[max_end_i]);
	  if (max_end_i > 0) value *= (channel == 0 ? gamma[max_end_i] : gamma[max_end_i] + beta[max_end_i] - gamma[max_end_i] * beta[max_end_i]);
	  psum = std::max(psum - value, 0.0);
	}
	psum = (channel == 0 ? psum * (1.0 - gamma[pos + 1]) + gamma[pos + 1]: psum * (1.0 - gamma[pos + 1]) * (1.0 - beta[pos + 1]) + (gamma[pos + 1] + beta[pos + 1] - gamma[pos + 1] * beta[pos + 1]));
      }
    }
    
    // M step
    double dc, cc; // dc: drop-off count; cc: covering count
    
    start2[0] += start[0];
    end2[0] += end[0];
    for (int i = 1; i <= len; ++i) {
      start2[i] += start[i] + start2[i - 1];
      end2[i] += end[i] + end2[i - 1];
      
      dc = std::max(0.0, end2[i] - end2[i - 1]); // drop-off count
      cc = std::max(0.0, end2[i] - start2[i - 1] - dc); // covering cout
      
      switch(getState()) {
      case 0: 
	// learn separately, (-) channel 
	if (isMAP) {
	  gamma[i] = (dgamma + dc) / (dgamma + dc + cgamma + cc);
	}
	else {
	  gamma[i] = (dc > 0.0 ? dc / (dc + cc) : 0.0);
	}
	break;
      case 1:
	// learn separately, (+) channel
	if (isMAP) {
	  solveQuadratic1(beta[i], gamma[i], dc, cc);
	}
	else {
	  beta[i] = (dc > 0.0 ? dc / (dc + cc) : 0.0);
	  beta[i] = ((beta[i] > gamma[i]) && (gamma[i] < 1.0) ? (beta[i] - gamma[i]) / (1.0 - gamma[i]) : 0.0);
	  if (isZero(1.0 - beta[i])) beta[i] = 1.0 - 1e-8; // truncate beta to be < 1 to calculate crate
	}
	break;
      case 2:
	dcm[i] = dc;
	ccm[i] = cc;
	break;
      case 3:
	if (isMAP) {
	  solveQuadratic2(gamma[i], beta[i], dcm[i], ccm[i], dc, cc);
	}
	else {
	  gamma[i] = (dcm[i] > 0.0 ? dcm[i] / (dcm[i] + ccm[i]) : 0.0);
	  beta[i] = (dc > 0.0 ? dc / (dc + cc) : 0.0);
	  if (beta[i] > gamma[i]) { 
	    assert(gamma[i] < 1.0);
	    beta[i] = (beta[i] - gamma[i]) / (1.0 - gamma[i]);
	  }
	  else {
	    gamma[i] = (dcm[i] + dc > 0.0 ? (dcm[i] + dc) / (dcm[i] + ccm[i] + dc + cc) : 0.0);
	    beta[i] = 0.0;
	  }
	  if (isZero(1.0 - beta[i])) beta[i] = 1.0 - 1e-8; // truncate beta to be < 1 to calculate crate
	}
	break;
      default: assert(false);
      }
    }
  }

  // Prepare for the next round
  calcAuxiliaryArrays(isJoint()? channel ^ 1 : channel);
}

void DMSTransModel::read(std::ifstream& fin, int channel) {
  std::string tmp_name;
  int tmp_len;

  fin>> tmp_name>> tmp_len;

  if (name == "") { 
    name = tmp_name; len = tmp_len;
    efflen = len - min_frag_len + 1;
    delta = 1.0 / (len + 1);
    if (efflen > 0) {
      // auxiliary arrays
      logsum = new double[len + 1];
      margin_prob = new double[efflen];
    }
  }
  else assert((tmp_name == name) && (tmp_len == len));

  if (channel == 0) {
    if (gamma == NULL) gamma = new double[len + 1];
    gamma[0] = 0.0;
    for (int i = 1; i <= len; ++i) fin>> gamma[i];
  }
  else {
    if (beta == NULL) beta = new double[len + 1];
    beta[0] = 0.0;
    for (int i = 1; i <= len; ++i) fin>> beta[i];
  }
}

void DMSTransModel::write(std::ofstream& fout, int channel) {
  fout<< name<< '\t'<< len;

  if (channel == 0) {
    for (int i = 1; i <= len; ++i) fout<< '\t'<< gamma[i];
  }
  else {
    for (int i = 1; i <= len; ++i) fout<< '\t'<< beta[i];
  }

  fout<< std::endl;
}

void DMSTransModel::writeFreq(std::ofstream& fc, std::ofstream& fout) {
  double c = 0.0;
  for (int i = 1; i <= len; ++i) c += -log(1.0 - beta[i]);

  fc << c;
  fout<< name;
  if (!isZero(c)) {
    fout<< '\t'<< len;
    for (int i = 1; i <= len; ++i) fout<< '\t'<< std::max(0.0, -log(1.0 - beta[i]) / c);
  }
  fout<< std::endl;
}

void DMSTransModel::startSimulation() {
  if (efflen <= 0) return;
  cdf_end = new double[efflen];

  calcAuxiliaryArrays(getChannel());
 
  cdf_end[0] = exp(logsum[min_frag_len] - logsum[0]) * margin_prob[0];
  for (int i = 1; i < efflen; ++i) {
    cdf_end[i] = (getChannel() == 0 ? gamma[i] : (gamma[i] + beta[i] - gamma[i] * beta[i])) * exp(logsum[i + min_frag_len] - logsum[i]) * margin_prob[i];
    cdf_end[i] += cdf_end[i - 1];
  }
}

void DMSTransModel::simulate(Sampler* sampler, int& pos, int& fragment_length) {
  int upper_bound, cpos;
  double random_value, value, sum;

  // Determine the end position
  pos = sampler->sample(cdf_end, efflen);

  // Determine fragment length
  upper_bound = std::min(max_frag_len, len - pos);
  random_value = sampler->random() * margin_prob[pos];
  value = sum = 1.0; 
  cpos = pos + min_frag_len;
  // [ ) intervals
  for (fragment_length = min_frag_len; (fragment_length < upper_bound) && (sum <= random_value); ++fragment_length) {
    ++cpos;
    value *= (getChannel() == 0 ? (1.0 - gamma[cpos]) : (1.0 - gamma[cpos]) * (1.0 - beta[cpos])); 
    sum += value;
  }
  
  fragment_length += primer_length;
}

void DMSTransModel::finishSimulation() {
  if (cdf_end != NULL) delete[] cdf_end;
  cdf_end = NULL;
}

// private member function
void DMSTransModel::calcAuxiliaryArrays(int channel) {
  double value;
  int max_pos;

  // Calculate logsum
  logsum[0] = 0.0;
  for (int i = 1; i <= len; ++i) {
    if (gamma[i] >= 1.0 || (channel == 1 && beta[i] >= 1.0)) value = -INF;
    else value = (channel == 0 ? log(1.0 - gamma[i]) : log(1.0 - gamma[i]) + log(1.0 - beta[i]));
    logsum[i] = logsum[i - 1] + value;
  }

  // Calculate margin_prob
  margin_prob[efflen - 1] = 1.0;
  for (int i = efflen - 2, pos = len; i >= 0; --i, --pos) {
    max_pos = (i + 1) + max_frag_len;
    assert(max_pos > len || margin_prob[i + 1] - exp(logsum[max_pos] - logsum[pos]) >= 0.0);
    margin_prob[i] = 1.0 + (channel == 0 ? (1.0 - gamma[pos]) : (1.0 - gamma[pos]) * (1.0 - beta[pos])) * \
      (max_pos > len ? margin_prob[i + 1] : margin_prob[i + 1] - exp(logsum[max_pos] - logsum[pos]));
  }

  // Calculate the probability of passing the size selection step
  prob_pass[channel] = 0.0;
  for (int i = 0; i < efflen; ++i) {
    value = delta * margin_prob[i] * exp(logsum[i + min_frag_len] - logsum[i]);
    if (i > 0) value *= (channel == 0 ? gamma[i] : gamma[i] + beta[i] - gamma[i] * beta[i]);
    prob_pass[channel] += value;
  }

  if (efflen != efflen2 && efflen2 > 0) {
    // Calculate marginal probability array for allocating SE reads 
    margin_prob2[efflen2 - 1] = 1.0;
    for (int i = efflen2 - 2, pos = len; i >= 0; --i, --pos) {
      max_pos = i + max_frag_len + 1;
      assert(max_pos > len || margin_prob2[i + 1] - exp(logsum[max_pos] - logsum[pos]) >= 0.0);
      margin_prob2[i] = 1.0 + (channel == 0 ? (1.0 - gamma[pos]) : (1.0 - gamma[pos]) * (1.0 - beta[pos])) * \
	(max_pos > len ? margin_prob2[i + 1] : margin_prob2[i + 1] - exp(logsum[max_pos] - logsum[pos]));
    }
  }
}
