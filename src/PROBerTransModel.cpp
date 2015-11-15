#include<cmath>
#include<cstring>
#include<cassert>
#include<string>
#include<fstream>
#include<algorithm>

#include "utils.h"
#include "sampling.hpp"
#include "PROBerTransModel.hpp"

const double PROBerTransModel::INF = 1000.0;

int PROBerTransModel::primer_length = 6; // default, 6bp
int PROBerTransModel::min_frag_len;
int PROBerTransModel::max_frag_len;

double PROBerTransModel::gamma_init;
double PROBerTransModel::beta_init;

double PROBerTransModel::base;
double PROBerTransModel::cgamma;
double PROBerTransModel::dgamma;
double PROBerTransModel::cbeta;
double PROBerTransModel::dbeta;

double PROBerTransModel::lgammas[2];
double PROBerTransModel::defaults[2];

int PROBerTransModel::min_alloc_len;
bool PROBerTransModel::isMAP = true; // default is true

bool PROBerTransModel::learning = false; // default is simulation
int PROBerTransModel::state = 0;

double PROBerTransModel::prob_p = 1.0; // default is no enrichment for signal
bool PROBerTransModel::enrich4signal = false; // default is not to enrich for signal
  
void PROBerTransModel::setGlobalParams(int primer_length, int min_frag_len, int max_frag_len, int init_state) { 
  assert(primer_length <= min_frag_len && min_frag_len <= max_frag_len);
  PROBerTransModel::primer_length = primer_length;
  PROBerTransModel::min_frag_len = min_frag_len - primer_length;
  PROBerTransModel::max_frag_len = max_frag_len - primer_length;
  state = init_state;
}

void PROBerTransModel::setLearningRelatedParams(double gamma_init, double beta_init, double base, int read_length, bool isMAP, bool enrich4signal) {
  learning = true;

  PROBerTransModel::gamma_init = gamma_init;
  PROBerTransModel::beta_init = beta_init;
  PROBerTransModel::base = base;
  PROBerTransModel::min_alloc_len = std::max(min_frag_len, read_length - primer_length);
  PROBerTransModel::isMAP = isMAP;
  PROBerTransModel::enrich4signal = enrich4signal;
  
  if (isMAP) {
    dgamma = gamma_init * base;
    cgamma = base - dgamma;
    dbeta = beta_init * base;
    cbeta = base - dbeta;

    lgammas[0] = lgamma(base + 2.0) - lgamma(dgamma + 1.0) - lgamma(cgamma + 1.0);
    lgammas[1] = lgamma(base + 2.0) - lgamma(dbeta + 1.0) - lgamma(cbeta + 1.0);
    defaults[0] = dgamma * log(gamma_init) + cgamma * log(1.0 - gamma_init);
    defaults[1] = dbeta * log(beta_init) + cbeta * log(1.0 - beta_init);
  }
}



PROBerTransModel::PROBerTransModel(int tid, const std::string& name, int transcript_length) : tid(tid), name(name) {  
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

  log_prior[0] = log_prior[1] = 0.0;

  cdf_end = NULL;

  if (!learning) return;

  len = transcript_length - primer_length;
  efflen = len - min_frag_len + 1;
  delta = 1.0 / (len + 1); // either len + 1 for both random priming and fragmentation or len for both; for now, len + 1, later will try len
  
  gamma = new double[len + 1];
  gamma[0] = 0.0;
  for (int i = 1; i <= len; ++i) gamma[i] = (isMAP ? gamma_init : 0.0);

  if (getState() > 0) {
    beta = new double[len + 1];
    beta[0] = 0.0;
    for (int i = 1; i <= len; ++i) beta[i] = (isMAP ? beta_init : 0.0);
  }
}

PROBerTransModel::~PROBerTransModel() {
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

void PROBerTransModel::init() {
  int state = getState();

  assert(state < 3);

  // set initial values for EM
  if (state != 1) for (int i = 1; i <= len; ++i) gamma[i] = gamma_init;
  if (state != 0) for (int i = 1; i <= len; ++i) beta[i] = beta_init;

  // count vectors for fragment starts and ends
  start = new double[len + 1];
  end = new double[len + 1];
  memset(start, 0, sizeof(double) * (len + 1));
  memset(end, 0, sizeof(double) * (len + 1));

  // Auxiliary arrays
  logsum = new double[len + 1];
  margin_prob = new double[efflen];
  memset(logsum, 0, sizeof(double) * (len + 1));
  memset(margin_prob, 0, sizeof(double) * efflen);
  
  if (state == 2) {
    dcm = new double[len + 1];
    ccm = new double[len + 1];
    memset(dcm, 0, sizeof(double) * (len + 1));
    memset(ccm, 0, sizeof(double) * (len + 1));
  }
  
  if (hasSE) {
    end_se = new double[len + 1];
    memset(end_se, 0, sizeof(double) * (len + 1));

    efflen2 = len - min_alloc_len + 1;
    assert(efflen2 > 0);
    if (efflen2 == efflen) efflen2 = -1; // If equal, do not need to build margin_prob2

    if (efflen2 > 0) {
      margin_prob2 = new double[efflen2];
      memset(margin_prob2, 0, sizeof(double) * (efflen2));
    }
  }
}

void PROBerTransModel::calcAuxiliaryArrays(int channel) {
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

  // Calculate the probability of passing the size selection step and the enrichment step
  prob_pass[channel] = 0.0;
  for (int i = 0; i < efflen; ++i) {
    value = delta * margin_prob[i] * exp(logsum[i + min_frag_len] - logsum[i]);
    if (channel == 0) { if (i > 0) value *= gamma[i]; }
    else { value *= (i > 0 ? beta[i] + (1.0 - beta[i]) * gamma[i] * prob_p : prob_p); } 
    prob_pass[channel] += value;
  }

  if (efflen2 > 0) {
    // Calculate marginal probability array for allocating SE reads 
    margin_prob2[efflen2 - 1] = 1.0;
    for (int i = efflen2 - 2, pos = len; i >= 0; --i, --pos) {
      max_pos = i + max_frag_len + 1;
      assert(max_pos > len || margin_prob2[i + 1] - exp(logsum[max_pos] - logsum[pos]) >= 0.0);
      margin_prob2[i] = 1.0 + (channel == 0 ? (1.0 - gamma[pos]) : (1.0 - gamma[pos]) * (1.0 - beta[pos])) * \
	(max_pos > len ? margin_prob2[i + 1] : margin_prob2[i + 1] - exp(logsum[max_pos] - logsum[pos]));
    }
  }

  if (isMAP) {
    log_prior[channel] = 0.0;
    if (channel == 0)
      for (int i = 1; i <= len; ++i) log_prior[channel] += dgamma * log(gamma[i]) + cgamma * log (1.0 - gamma[i]);
    else 
      for (int i = 1; i <= len; ++i) log_prior[channel] += dbeta * log(beta[i]) + cbeta * log(1.0 - beta[i]);
  }
}

void PROBerTransModel::update() {
  int channel = getChannel();
  std::vector<InMemAlign*> &alignments = alignmentsArr[channel];

  HIT_INT_TYPE size = alignments.size();

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

void PROBerTransModel::EM_step(double N_tot, double& c_4_p, double& c_4_1mp) {
  int channel = getChannel();

  int max_end_i;
  double psum, value, ecpp; // ecpp, expected count per position

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

    ecpp = N_tot * delta;

    //E step, if we have reads that do not know their start positions, infer start from end
    if (!isZero(N_se)) {
      double curr;
      double *mp = NULL;
      int effl;

      if (min_alloc_len > min_frag_len) { mp = margin_prob2; effl = efflen2; }
      else { mp = margin_prob; effl = efflen; }
      assert(effl > 0);

      curr = (end_se[0] > 0.0 && mp[0] > 0.0) ? end_se[0] / mp[0] : 0.0;
      start[min_alloc_len] += curr;
      for (int i = 1, pos = min_alloc_len + 1; i < effl; ++i, ++pos) {
	value = curr;
	curr = (end_se[i] > 0.0 && mp[i] > 0.0) ? end_se[i] / mp[i] : 0.0;

	max_end_i = (pos - 1) - max_frag_len;	
	if (max_end_i >= 0) {
	  value -= ((end_se[max_end_i] > 0.0 && mp[max_end_i] > 0.0) ? end_se[max_end_i] * (exp(logsum[pos - 1] - logsum[max_end_i + min_alloc_len]) / mp[max_end_i]) : 0.0);
	  if (value < 0.0) value = 0.0;
	}
	
	curr += (channel == 0 ? (1.0 - gamma[pos]) : (1.0 - gamma[pos]) * (1.0 - beta[pos])) * value;
	start[pos] += curr;
      }
    }



    // E step, calculate hidden reads
    // Calculate end2, do not count for reads that fail to enrich
    psum = 1.0; // sum of probability of starting from anywhere and end at i, with the drop-off probability excluded
    for (int i = len; i >= 0; --i) {
      end2[i] = psum;
      if (i < efflen) end2[i] = std::max(end2[i] - exp(logsum[i + min_frag_len] - logsum[i]) * margin_prob[i], 0.0);
      value = (i > 0 ? (channel == 0 ? gamma[i] : beta[i] + (1.0 - beta[i]) * gamma[i]) : 1.0); // drop-off probability for psum
      end2[i] *= ecpp * value;
      if (i > 0) psum = 1.0 + psum * (channel == 0 ? (1.0 - gamma[i]) : (1.0 - gamma[i]) * (1.0 - beta[i]));
    }


    
    // Calculate start2
    for (int i = 0; i < min_frag_len; ++i) start2[i] = ecpp;
    psum = (channel == 0 ? 1.0 : prob_p); // sum of size-selection-passed and enriched fragments, the min_frag_len portion excluded
    for (int i = min_frag_len, pos = 0; i <= len; ++i, ++pos) {
      start2[i] = std::max(1.0 - psum * exp(logsum[i] - logsum[pos]), 0.0);
      start2[i] *= ecpp;
      if (i < len) {
	max_end_i = i - max_frag_len;
	if (max_end_i >= 0) {
	  value = exp(logsum[pos] - logsum[max_end_i]);
	  if (channel == 0) { if (max_end_i > 0) value *= gamma[max_end_i]; }
	  else { value *= (max_end_i > 0 ? beta[max_end_i] + (1.0 - beta[max_end_i]) * gamma[max_end_i] * prob_p : prob_p); }
	  psum = std::max(psum - value, 0.0);
	}
	psum = (channel == 0 ? psum * (1.0 - gamma[pos + 1]) + gamma[pos + 1]: \
		psum * (1.0 - gamma[pos + 1]) * (1.0 - beta[pos + 1]) + (beta[pos + 1] + (1.0 - beta[pos + 1]) * gamma[pos + 1] * prob_p));
      }
    }


    
    // Collect sufficient statistics for prob_p
    if (enrich4signal && channel == 1) {
      double prob_1mp = 1.0 - prob_p;
      
      c_4_p += end[0];
      c_4_1mp += ecpp * margin_prob[0] * exp(logsum[min_frag_len]) * prob_1mp;
      for (int i = 1; i < efflen; ++i) {
	value = (1.0 - beta[i]) * gamma[i];
	if (end[i] > 0.0) c_4_p += end[i] * (value * prob_p / (beta[i] + value * prob_p));
	c_4_1mp += ecpp * margin_prob[i] * exp(logsum[i + min_frag_len] - logsum[i]) * value * prob_1mp;
      }
    }


    
    // M step
    double dc, cc, d1; // dc: drop-off count; cc: covering count; d1: drop-of count due to modification
    
    start2[0] += start[0];
    end2[0] += end[0];
    for (int i = 1; i <= len; ++i) {

      // calculate D1
      if (channel == 1) {
	d1 = 0.0;
	if (end[i] > 0.0) d1 += end[i] * (beta[i] / (beta[i] + (1.0 - beta[i]) * gamma[i] * prob_p));
	if (end2[i] > 0.0) d1 += end2[i] * (beta[i] / (beta[i] + (1.0 - beta[i]) * gamma[i]));
	if (i < efflen) end2[i] += ecpp * margin_prob[i] * exp(logsum[i + min_frag_len] - logsum[i]) * (1.0 - beta[i]) * gamma[i] * (1.0 - prob_p);
      }
      
      end2[i] += end[i];
      dc = end2[i]; // drop-off count

      start2[i] += start[i] + start2[i - 1];
      end2[i] += end2[i - 1];

      cc = std::max(0.0, end2[i] - start2[i - 1] - dc); // covering cout
      
      switch(getState()) {
      case 0: 
	// learn separately, (-) channel 
	if (isMAP) {
	  gamma[i] = (dgamma + dc) / (dgamma + dc + cgamma + cc);
	  assert(gamma[i] > 0.0 && gamma[i] < 1.0);
	}
	else
	  gamma[i] = (dc > 0.0 ? dc / (dc + cc) : 0.0);
	
	break;
      case 1:
	// learn separately, (+) channel
	assert(dc - d1 >= 0.0);
	
	if (isMAP) {
	  beta[i] = (dbeta + d1) / (dbeta + dc + cbeta + cc);
	  assert(beta[i] > 0.0 && beta[i] < 1.0);
	}
	else 
	  beta[i] = (d1 > 0.0 ? d1 / (dc + cc) : 0.0);
	
	break;
      case 2:
	dcm[i] = dc;
	ccm[i] = cc;
	
	break;
      case 3:
	assert(dc - d1 >= 0.0);
	
	if (isMAP) {
	  gamma[i] = (dgamma + dcm[i] + (dc - d1)) / (dgamma + dcm[i] + (dc - d1) + cgamma + ccm[i] + cc);
	  beta[i] = (dbeta + d1) / (dbeta + dc + cbeta + cc);
	  assert(gamma[i] > 0.0 && gamma[i] < 1.0 && beta[i] > 0.0 && beta[i] < 1.0);
	}
	else {
	  gamma[i] = (dcm[i] + (dc - d1) > 0.0 ? (dcm[i] + (dc - d1)) / (dcm[i] + (dc - d1) + ccm[i] + cc) : 0.0);
	  beta[i] = (d1 > 0.0 ? d1 / (dc + cc) : 0.0);
	}
	
	break;
      default: assert(false);
      }
    }
  }
  
  // Prepare for the next round
  calcAuxiliaryArrays(isJoint()? channel ^ 1 : channel);
}

void PROBerTransModel::read(std::ifstream& fin, int channel) {
  std::string tmp_name;
  int tmp_len;

  fin>> tmp_name>> tmp_len;

  if (name == "") { 
    name = tmp_name; len = tmp_len;
    efflen = len - min_frag_len + 1;
    delta = 1.0 / (len + 1); // either len + 1 for both random priming and fragmentation or len for both; for now, len + 1, later will try len
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

void PROBerTransModel::write(std::ofstream& fout, int channel) {
  fout<< name<< '\t'<< len;

  if (channel == 0) {
    for (int i = 1; i <= len; ++i) fout<< '\t'<< gamma[i];
  }
  else {
    for (int i = 1; i <= len; ++i) fout<< '\t'<< beta[i];
  }

  fout<< std::endl;
}

void PROBerTransModel::startSimulation() {
  if (efflen <= 0) return;
  cdf_end = new double[efflen];

  calcAuxiliaryArrays(getChannel());
 
  cdf_end[0] = prob_p * exp(logsum[min_frag_len]) * margin_prob[0];
  for (int i = 1; i < efflen; ++i) {
    cdf_end[i] = (getChannel() == 0 ? gamma[i] : beta[i] + (1.0 - beta[i]) * gamma[i] * prob_p) * exp(logsum[i + min_frag_len] - logsum[i]) * margin_prob[i];
    cdf_end[i] += cdf_end[i - 1];
  }
}

void PROBerTransModel::simulate(Sampler* sampler, int& pos, int& fragment_length) {
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

void PROBerTransModel::finishSimulation() {
  if (cdf_end != NULL) delete[] cdf_end;
  cdf_end = NULL;
}
