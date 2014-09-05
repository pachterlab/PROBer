#include<cmath>
#include<cstring>
#include<cassert>
#include<string>
#include<fstream>
#include<algorithm>

#include "sampling.hpp"
#include "DMSTransModel.hpp"

DMSTransModel::DMSTransModel(bool learning, int transcript_length, Sampler* sampler) {
  gamma = beta = NULL;
  start = end = NULL;
  logsum = margin_prob = NULL;
  start2 = end2 = NULL;
  isSE = false;

  this->learning = learning;

  if (!learning) return;

  len = transcript_length - primer_length;
  efflen = len - min_frag_len + 1;
  delta = 1.0 / (len + 1.0);
  prob_pass = 0.0;
  
  gamma = new double[len + 1];
  memset(gamma, 0, sizeof(double) * (len + 1));

  if (efflen <= 0) return;

  // auxiliary arrays
  logsum = new double[len + 1];
  margin_prob = new double[efflen];

  for (int i = 1; i <= len; ++i) gamma[i] = sampler->sample();
  
  start = new double[len + 1];
  end = new double[len + 1];

  // Initialize arrays
  memset(start, 0, sizeof(double) * (len + 1));
  memset(end, 0, sizeof(double) * (len + 1));
  
  start2 = new double[len + 1];
  end2 = new double[len + 1];
}

DMSTransModel::DMSTransModel(const DMSTransModel& o) : learning(o.learning), isSE(o.isSE), len(o.len), efflen(o.efflen), delta(o.delta), prob_pass(o.prob_pass) {
  gamma = beta = NULL;
  start = end = NULL;
  logsum = margin_prob = NULL;
  start2 = end2 = NULL;

  if (o.gamma != NULL) {
    gamma = new double[len + 1];
    memcpy(gamma, o.gamma, sizeof(double) * (len + 1));
  }
  if (o.beta != NULL) {
    beta = new double[len + 1];
    memcpy(beta, o.beta, sizeof(double) * (len + 1));
  }
  if (o.start != NULL) {
    start = new double[len + 1];
    memcpy(start, o.start, sizeof(double) * (len + 1));
  }
  if (o.end != NULL) {
    end = new double[len + 1];
    memcpy(end, o.end, sizeof(double) * (len + 1));
  }
  if (o.logsum != NULL) {
    logsum = new double[len + 1];
    memcpy(logsum, o.logsum, sizeof(double) * (len + 1));
  }
  if (o.margin_prob != NULL) {
    margin_prob = new double[efflen];
    memcpy(margin_prob, o.margin_prob, sizeof(double) * efflen);
  }
  if (o.start2 != NULL) {
    start2 = new double[len + 1];
    memcpy(start2, o.start2, sizeof(double) * (len + 1));
  }
  if (o.end2 != NULL) {
    end2 = new double[len + 1];
    memcpy(end2, o.end2, sizeof(double) * (len + 1));
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

  // If shared across a thread of transcripts, no need to delete
  if (start2 != NULL) delete[] start2;
  if (end2 != NULL) delete[] end2;
}

const double DMSTransModel::eps = 1e-8;
const double DMSTransModel::INF = 1000.0;

int DMSTransModel::primer_length = 6; // default, 6bp
int DMSTransModel::min_frag_len;
int DMSTransModel::max_frag_len;

void DMSTransModel::setGlobalParams(int primer_length, int min_frag_len, int max_frag_len) {
  assert(primer_length <= min_frag_len && min_frag_len <= max_frag_len);
  DMSTransModel::primer_length = primer_length;
  DMSTransModel::min_frag_len = min_frag_len - primer_length;
  DMSTransModel::max_frag_len = max_frag_len - primer_length;
}

void DMSTransModel::calcAuxiliaryArrays() {
  double value;

  // Calculate logsum
  logsum[0] = 0.0;
  for (int i = 1; i <= len; ++i) {
    value = (beta == NULL ? log(1.0 - gamma[i]) : log((1.0 - gamma[i]) * (1.0 - beta[i])));
    if (isinf(value)) value = -INF;
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

void DMSTransModel::init() {
  memset(start, 0, sizeof(double) * (len + 1));
  memset(end, 0, sizeof(double) * (len + 1));
}

// May consider implement Kahan summation algorithm if precision is really a concern
void DMSTransModel::EM(double N_obs, int round) {
  int max_end_i;
  double psum, value;
  double N_tot = N_obs / prob_pass;

  for (int ROUND = 0; ROUND < round; ++ROUND) {
    // E step
    
    // Step E1, infer start from end if needed
    if (isSE) {
      start[min_frag_len] = end[0] / margin_prob[0];
      for (int i = 1, pos = min_frag_len + 1; i < efflen; ++i, ++pos) {
	start[pos] = end[i] / margin_prob[i];
	max_end_i = pos - max_frag_len - 1;
	start[pos] += (beta == NULL ? (1.0 - gamma[pos]) : (1.0 - gamma[pos]) * (1.0 - beta[pos])) * \
	  (max_end_i < 0 ? start[pos - 1] : start[pos - 1] - end[max_end_i] * (exp(logsum[pos - 1] - logsum[max_end_i + min_frag_len]) / margin_prob[max_end_i]));
	if (start[pos] <= 0.0) start[pos] = 0.0;
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
      value = std::max((end2[i] - end2[i - 1]) / (end2[i] - start2[i - 1]), 0.0);
      if (beta == NULL) gamma[i] = value;
      else {
	if (value - gamma[i] <= eps || 1.0 - gamma[i] <= eps) beta[i] = 0.0;
	else beta[i] = (value - gamma[i]) / (1.0 - gamma[i]);
      }
    }
    
    // Prepare for the next round
    calcAuxiliaryArrays();
    N_tot = N_obs / prob_pass;
  }
}

void DMSTransModel::read(std::ifstream& fin, Sampler* sampler) {
  int tmp_len;

  fin>> tmp_len;

  if (learning) {
    assert(tmp_len == len);
    if (beta == NULL) {
      gamma[0] = 0.0;
      for (int i = 1; i <= len; ++i) fin>> gamma[i];
      // Set initial values
      beta = new double[len + 1];
      memset(beta, 0, sizeof(double) * (len + 1));
      if (efflen > 0) for (int i = 1; i <= len; ++i) beta[i] = sampler->sample();
    }
    else {
      beta[0] = 0.0;
      for (int i = 1; i <= len; ++i) fin>> beta[i];
    }
  }
  else {
    if (gamma == NULL) {
      len = tmp_len;
      efflen = len - min_frag_len + 1;
      delta = 1.0 / (len + 1);
      prob_pass = 0.0;
      
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
      assert(beta == NULL && tmp_len == len);
      beta = new double[len + 1];
      beta[0] = 0.0;
      for (int i = 1; i <= len; ++i) fin>> beta[i];
    }
  }
}

void DMSTransModel::write(std::ofstream& fout) {
  fout<< len;

  fout.precision(10);
  fout.unsetf(std::ios::floatfield);

  if (beta == NULL) {
    for (int i = 1; i <= len; ++i) fout<< '\t'<< gamma[i];
  }
  else {
    for (int i = 1; i <= len; ++i) fout<< '\t'<< beta[i];
  }

  fout<< std::endl;
}

void DMSTransModel::writeTheta(std::ofstream& fout) {
  double c = 0.0;
  for (int i = 1; i <= len; ++i) c += -log(1.0 - beta[i]);

  fout.precision(10);
  fout.unsetf(std::ios::floatfield);

  fout<< c<< '\t'<< len;
  for (int i = 1; i <= len; ++i) fout<< '\t'<< std::max(0.0, -log(1.0 - beta[i]) / c);
  fout<< std::endl;
}

void DMSTransModel::simulate(Sampler* sampler, int& pos, int& fragment_length) {
  double value;

  do {
    pos = int(sampler->random() * (len + 1));
    fragment_length = 0;
    while (pos > 0) {
      value = sampler->random();
      if ((beta == NULL && value < gamma[pos]) || (beta != NULL && value < gamma[pos] + beta[pos] - gamma[pos] * beta[pos])) break;
      --pos;
      ++fragment_length;
    }
  } while (fragment_length < min_frag_len || fragment_length > max_frag_len);
  fragment_length += primer_length;
}
