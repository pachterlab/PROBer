#include<cmath>
#include<cstring>
#include<cassert>
#include<string>
#include<fstream>
#include<algorithm>

#include "utils.h"
#include "sampling.hpp"
#include "DMSTransModel.hpp"

DMSTransModel::DMSTransModel(bool learning, const std::string& name, int transcript_length) {  
  gamma = beta = NULL;
  start = end = NULL;

  logsum = margin_prob = NULL;
  margin_prob2 = NULL;

  start2 = end2 = NULL;
  alignments.clear();

  this->name = "";
  len = efflen = -1; 
  efflen2 = -1;
  N_obs = 0.0;
  delta = prob_pass = 0.0;

  cdf_end = NULL;

  this->learning = learning;

  if (!learning) return;

  this->name = name;
  len = transcript_length - primer_length;
  efflen = len - min_frag_len + 1;
  efflen2 = len - min_alloc_len + 1;
  delta = 1.0 / (len + 1.0);
  
  gamma = new double[len + 1];
  memset(gamma, 0, sizeof(double) * (len + 1));

  if (efflen <= 0) return;

  // auxiliary arrays
  logsum = new double[len + 1];
  margin_prob = new double[efflen];
  for (int i = 1; i <= len; ++i) gamma[i] = gamma_init;
  
  start = new double[len + 1];
  end = new double[len + 1];

  // Initialize arrays
  memset(start, 0, sizeof(double) * (len + 1));
  memset(end, 0, sizeof(double) * (len + 1));

  if (efflen != efflen2 && efflen2 > 0) margin_prob2 = new double[efflen2];
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
}

const double DMSTransModel::INF = 1000.0;

int DMSTransModel::primer_length = 6; // default, 6bp
int DMSTransModel::min_frag_len;
int DMSTransModel::max_frag_len;
double DMSTransModel::gamma_init;
double DMSTransModel::beta_init;
int DMSTransModel::min_alloc_len;

void DMSTransModel::setGlobalParams(int primer_length, int min_frag_len, int max_frag_len, double gamma_init, double beta_init, int read_length) {
  assert(primer_length <= min_frag_len && min_frag_len <= max_frag_len);
  DMSTransModel::primer_length = primer_length;
  DMSTransModel::min_frag_len = min_frag_len - primer_length;
  DMSTransModel::max_frag_len = max_frag_len - primer_length;

  DMSTransModel::gamma_init = gamma_init;
  DMSTransModel::beta_init = beta_init;

  DMSTransModel::min_alloc_len = (read_length < min_frag_len ? min_frag_len : read_length) - primer_length;
}

void DMSTransModel::calcAuxiliaryArrays() {
  double value;
  int max_pos;

  assert(efflen > 0);

  // Calculate logsum
  logsum[0] = 0.0;
  for (int i = 1; i <= len; ++i) {
    if (gamma[i] >= 1.0 || (beta != NULL && beta[i] >= 1.0)) value = -INF;
    else value = (beta == NULL ? log(1.0 - gamma[i]) : log((1.0 - gamma[i]) * (1.0 - beta[i])));
    logsum[i] = logsum[i - 1] + value;
  }

  // Calculate margin_prob
  margin_prob[efflen - 1] = 1.0;
  for (int i = efflen - 2, pos = len; i >= 0; --i, --pos) {
    max_pos = i + max_frag_len + 1;
    assert(max_pos > len || margin_prob[i + 1] - exp(logsum[max_pos] - logsum[pos]) >= 0.0);
    margin_prob[i] = 1.0 + (beta == NULL ? (1.0 - gamma[pos]) : (1.0 - gamma[pos]) * (1.0 - beta[pos])) * \
      (max_pos > len ? margin_prob[i + 1] : margin_prob[i + 1] - exp(logsum[max_pos] - logsum[pos]));
  }

  // Calculate the probability of passing the size selection step
  prob_pass = 0.0;
  for (int i = 0; i < efflen; ++i) {
    value = delta * margin_prob[i] * exp(logsum[i + min_frag_len] - logsum[i]);
    if (i > 0) value *= (beta == NULL ? gamma[i] : gamma[i] + beta[i] - gamma[i] * beta[i]);
    prob_pass += value;
  }

  if (efflen != efflen2 && efflen2 > 0) {
    // Calculate marginal probability array for allocating SE reads 
    margin_prob2[efflen2 - 1] = 1.0;
    for (int i = efflen2 - 2, pos = len; i >= 0; --i, --pos) {
      max_pos = i + max_frag_len + 1;
      assert(max_pos > len || margin_prob2[i + 1] - exp(logsum[max_pos] - logsum[pos]) >= 0.0);
      margin_prob2[i] = 1.0 + (beta == NULL ? (1.0 - gamma[pos]) : (1.0 - gamma[pos]) * (1.0 - beta[pos])) * \
	(max_pos > len ? margin_prob2[i + 1] : margin_prob2[i + 1] - exp(logsum[max_pos] - logsum[pos]));
    }
  }
}

void DMSTransModel::update() {
  HIT_INT_TYPE size = alignments.size();
  double N_se;

  if (size == 0) return;

  // initialize
  N_obs = 0.0;
  memset(start, 0, sizeof(double) * (len + 1));
  memset(end, 0, sizeof(double) * (len + 1));

  assert(start2 != NULL && end2 != NULL);
  memset(end2, 0, sizeof(double) * (len + 1));

  N_se = 0.0;
  for (HIT_INT_TYPE i = 0; i < size; ++i) {
    end[alignments[i]->pos] += alignments[i]->frac;
    if (alignments[i]->fragment_length > 0) {
      start[alignments[i]->pos + alignments[i]->fragment_length - primer_length] += alignments[i]->frac; 
    }
    else {
      end2[alignments[i]->pos] += alignments[i]->frac;
      N_se += alignments[i]->frac;
    }
    N_obs += alignments[i]->frac;
  }    

  if (!isZero(N_se)) {
    // Infer start from end
    int max_end_i;    
    double prev, curr, value;
    double *mp = NULL;

    assert(efflen2 > 0);

    mp = (min_alloc_len > min_frag_len ? margin_prob2 : margin_prob);
    curr = (end2[0] > 0.0 && mp[0] > 0.0) ? end2[0] / mp[0] : 0.0;
    start[min_alloc_len] += curr;
    for (int i = 1, pos = min_alloc_len + 1; i < efflen2; ++i, ++pos) {
      prev = curr;
      curr = (end2[i] > 0.0 && mp[i] > 0.0) ? end2[i] / mp[i] : 0.0;
      max_end_i = pos - max_frag_len - 1;

      value = prev;
      if (max_end_i >= 0) {
	value -= ((end2[max_end_i] > 0.0 && mp[max_end_i] > 0.0) ? end2[max_end_i] * (exp(logsum[pos - 1] - logsum[max_end_i + min_alloc_len]) / mp[max_end_i]) : 0.0);
	if (value < 0.0) value = 0.0;
      }

      curr += (beta == NULL ? (1.0 - gamma[pos]) : (1.0 - gamma[pos]) * (1.0 - beta[pos])) * value;
      start[pos] += curr;
    }
  }
}

void DMSTransModel::EM(double N_tot) {
  if (efflen <= 0 || isZero(N_obs)) return;

  int max_end_i;
  double psum, value;
  //  double N_tot = N_obs / prob_pass;

  assert(start2 != NULL && end2 != NULL);

  // E step, calculate hidden reads
  // Calculate end2
  psum = 1.0;
  for (int i = len; i >= 0; --i) {
    if (i < efflen) end2[i] = std::max(psum - exp(logsum[i + min_frag_len] - logsum[i]) * margin_prob[i], 0.0);
    else end2[i] = psum;
    if (i > 0) end2[i] *= (beta == NULL ? gamma[i] : gamma[i] + beta[i] - gamma[i] * beta[i]);
    end2[i] *= delta * N_tot;
    if (i > 0) psum = 1.0 + psum * (beta == NULL ? (1.0 - gamma[i]) : (1.0 - gamma[i]) * (1.0 - beta[i]));
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
	if (max_end_i > 0) value *= (beta == NULL ? gamma[max_end_i] : gamma[max_end_i] + beta[max_end_i] - gamma[max_end_i] * beta[max_end_i]);
	psum = std::max(psum - value, 0.0);
      }
      psum = (beta == NULL ? psum * (1.0 - gamma[pos + 1]) + gamma[pos + 1]: psum * (1.0 - gamma[pos + 1]) * (1.0 - beta[pos + 1]) + (gamma[pos + 1] + beta[pos + 1] - gamma[pos + 1] * beta[pos + 1]));
    }
  }
  
  // M step
  start2[0] += start[0];
  end2[0] += end[0];
  for (int i = 1; i <= len; ++i) {
    start2[i] += start[i] + start2[i - 1];
    end2[i] += end[i] + end2[i - 1];      
    value = (end2[i] > end2[i - 1]) && (end2[i] > start2[i - 1]) ? (end2[i] - end2[i - 1]) / (end2[i] - start2[i - 1]) : 0.0;
    if (beta == NULL) gamma[i] = isZero(1.0 - value) ? 1.0 - 1e-8 : value; // avoid insufficient data lead to failing rate of 1, may seek better estimator in the future
    else { 
      beta[i] = (value > gamma[i]) && (gamma[i] < 1.0) ? (value - gamma[i]) / (1.0 - gamma[i]) : 0.0; 
      if (isZero(1.0 - beta[i])) beta[i] = 1.0 - 1e-8; // if floating point inaccuracy leads to beta > 1.0, reset it to be 1.0 - 1e-8
    }
  }
  
  // Prepare for the next round
  calcAuxiliaryArrays();
}

void DMSTransModel::read(std::ifstream& fin) {
  std::string tmp_name;
  int tmp_len;

  fin>> tmp_name>> tmp_len;

  if (learning) {
    assert((tmp_name == name) && (tmp_len == len));
    if (beta == NULL) {
      gamma[0] = 0.0;
      for (int i = 1; i <= len; ++i) fin>> gamma[i];
      // Set initial values
      beta = new double[len + 1];
      memset(beta, 0, sizeof(double) * (len + 1));
      if (efflen > 0) for (int i = 1; i <= len; ++i) beta[i] = beta_init;
    }
    else {
      assert(false);
      beta[0] = 0.0;
      for (int i = 1; i <= len; ++i) fin>> beta[i];
    }
  }
  else {
    if (gamma == NULL) {
      name = tmp_name;
      len = tmp_len;
      efflen = len - min_frag_len + 1;
      delta = 1.0 / (len + 1);
      
      gamma = new double[len + 1];
      gamma[0] = 0.0;
      for (int i = 1; i <= len; ++i) fin>> gamma[i];

      if (efflen > 0) {
	// auxiliary arrays
	logsum = new double[len + 1];
	margin_prob = new double[efflen];
      } 
    }
    else {
      assert((beta == NULL) && (tmp_name == name) && (tmp_len == len));
      beta = new double[len + 1];
      beta[0] = 0.0;
      for (int i = 1; i <= len; ++i) fin>> beta[i];
    }
  }
}

void DMSTransModel::write(std::ofstream& fout) {
  fout<< name<< '\t'<< len;

  if (beta == NULL) {
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

  calcAuxiliaryArrays();
 
  cdf_end[0] = exp(logsum[min_frag_len] - logsum[0]) * margin_prob[0];
  for (int i = 1; i < efflen; ++i) {
    cdf_end[i] = (beta == NULL ? gamma[i] : (gamma[i] + beta[i] - gamma[i] * beta[i])) * exp(logsum[i + min_frag_len] - logsum[i]) * margin_prob[i];
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
  for (fragment_length = min_frag_len; (fragment_length < upper_bound) && (sum <= random_value); ++fragment_length) {
    ++cpos;
    value *= (beta == NULL ? (1.0 - gamma[cpos]) : (1.0 - gamma[cpos]) * (1.0 - beta[cpos])); 
    sum += value;
  }
  
  fragment_length += primer_length;
}

void DMSTransModel::finishSimulation() {
  if (cdf_end != NULL) delete[] cdf_end;
  cdf_end = NULL;
}
