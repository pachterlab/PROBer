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
  start2 = end2 = NULL;
  isSE = false;
  alignments.clear();

  this->name = "";
  len = efflen = -1; 
  N_obs = 0.0;
  delta = prob_pass = 0.0;
  
  cdf_end = NULL;

  this->learning = learning;

  if (!learning) return;

  this->name = name;
  len = transcript_length - primer_length;
  efflen = len - min_frag_len + 1;
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
}

const double DMSTransModel::INF = 1000.0;

int DMSTransModel::primer_length = 6; // default, 6bp
int DMSTransModel::min_frag_len;
int DMSTransModel::max_frag_len;
double DMSTransModel::gamma_init;
double DMSTransModel::beta_init;

void DMSTransModel::setGlobalParams(int primer_length, int min_frag_len, int max_frag_len, double gamma_init, double beta_init) {
  assert(primer_length <= min_frag_len && min_frag_len <= max_frag_len);
  DMSTransModel::primer_length = primer_length;
  DMSTransModel::min_frag_len = min_frag_len - primer_length;
  DMSTransModel::max_frag_len = max_frag_len - primer_length;

  DMSTransModel::gamma_init = gamma_init;
  DMSTransModel::beta_init = beta_init;
}

void DMSTransModel::calcAuxiliaryArrays() {
  double value;

  // Calculate logsum
  logsum[0] = 0.0;
  for (int i = 1; i <= len; ++i) {
    if (gamma[i] >= 1.0 || (beta != NULL && beta[i] >= 1.0)) value = -INF;
    else value = (beta == NULL ? log(1.0 - gamma[i]) : log((1.0 - gamma[i]) * (1.0 - beta[i])));
    logsum[i] = logsum[i - 1] + value;
  }

  // Calculate margin_prob
  margin_prob[efflen - 1] = 1.0;
  int max_pos;
  for (int i = efflen - 2, pos = len; i >= 0; --i, --pos) {
    max_pos = i + max_frag_len + 1;
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
}

void DMSTransModel::EM(double N_tot, int round) {
  if (efflen <= 0 || isZero(N_obs)) return;

  int max_end_i;
  double psum, value;
  //  double N_tot = N_obs / prob_pass;

  assert(round == 1);

  assert(start2 != NULL && end2 != NULL);

  for (int ROUND = 0; ROUND < round; ++ROUND) {
    // E step
    
    // Step E1, infer start from end if needed
    if (isSE) {
      start[min_frag_len] = (end[0] > 0.0 && margin_prob[0] > 0.0) ? end[0] / margin_prob[0] : 0.0;
      for (int i = 1, pos = min_frag_len + 1; i < efflen; ++i, ++pos) {
	start[pos] = (end[i] > 0.0 && margin_prob[i] > 0.0) ? end[i] / margin_prob[i] : 0.0;
	max_end_i = pos - max_frag_len - 1;
	start[pos] += (beta == NULL ? (1.0 - gamma[pos]) : (1.0 - gamma[pos]) * (1.0 - beta[pos])) * \
	  (max_end_i < 0 ? start[pos - 1] : start[pos - 1] - \
	   ((end[max_end_i] > 0.0 && margin_prob[max_end_i] > 0.0) ? end[max_end_i] * (exp(logsum[pos - 1] - logsum[max_end_i + min_frag_len]) / margin_prob[max_end_i]) : 0.0));
      }
    }

    // Step E2, calculate hidden reads
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
      else beta[i] = (value > gamma[i]) && (gamma[i] < 1.0) ? (value - gamma[i]) / (1.0 - gamma[i]) : 0.0;
    }
    
    // Prepare for the next round
    calcAuxiliaryArrays();
  }
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
