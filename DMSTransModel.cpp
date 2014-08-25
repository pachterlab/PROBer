#include<cstring>
#include<cassert>
#include<string>
#include<fstream>

#include "DMSTransModel.hpp"

DMSTransModel::DMSTransModel(int transcript_length, bool learning) {
  gamma = beta = NULL;
  start = end = NULL;
  logsum = margin_prob = NULL;
  start2 = end2 = NULL;
  isSE = false;

  len = transcript_length - primer_length;
  efflen = len - min_frag_len + 1;
  delta = 1.0 / (len + 1.0);
  prob_pass = 0.0;

  if (efflen <= 0) return;

  gamma = new double[len + 1];
  gamma[0] = 0.0;
  for (int i = 1; i <= len; ++i) gamma[i] = 0.01;

  // auxiliary arrays
  logsum = new double[len + 1];
  margin_prob = new double[efflen];

  if (!learning) return;
  
  start = new double[len + 1];
  end = new double[len + 1];

  // Initialize arrays
  memset(start, 0, sizeof(double) * (len + 1));
  memset(end, 0, sizeof(double) * (len + 1));
  
  start2 = new double[len + 1];
  end2 = new double[len + 1];
}

DMSTransModel::~DMSTransModel() {
  if (efflen <= 0) return;

  if (gamma != NULL) delete[] gamma;
  if (beta != NULL) delete[] beta;

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
    value = delta * margin_prob[i] *exp(logsum[i + min_frag_len] - logsum[i]);
    if (i > 0) value *= (beta == NULL ? gamma[i] : (gamma[i] + beta[i] - gamma[i] * beta[i]));
    prob_pass += value;
  }
}

void DMSTransModel::init() {
  memset(start, 0, sizeof(double) * (len + 1));
  memset(end, 0, sizeof(double) * (len + 1));
}

// May consider implement Kahan summation algorithm if precision is really a concern
void DMSTransModel::EM(double N_tot, int round) {
  int max_end_i;
  double psum, value;

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
	  psum -= value;
	}
	psum = (beta == NULL ? psum * (1.0 - gamma[pos + 1]) + gamma[pos + 1]: psum * (1.0 - gamma[pos + 1]) * (1.0 - beta[pos + 1]) + (gamma[pos] + beta[pos] - gamma[pos] * beta[pos]));
      }
    }

    // M step

    start2[0] += start[0];
    end2[0] += end[0];
    for (int i = 1; i <= len; ++i) {
      start2[i] += start[i] + start2[i - 1];
      end2[i] += end[i] + end2[i - 1];      
      value = (end2[i] - end2[i - 1]) / (end2[i] - start2[i - 1]);
      if (beta == NULL) gamma[i] = value;
      else {
	if (value - gamma[i] <= eps || gamma[i] - 1.0 <= eps) beta[i] = 0.0;
	else beta[i] = (value - gamma[i]) / (1.0 - gamma[i]);
      }
    }
    
    // Prepare for the next round
    calcAuxiliaryArrays();
  }
}

void DMSTransModel::read(std::ifstream fin) {
  int tmp_len;

  fin>> tmp_len;
  assert(tmp_len == len);
  
  if (beta == NULL) {
    gamma[0] = 0.0;
    for (int i = 1; i <= len; ++i) fin>> gamma[i];
    
    // Set initial values
    beta = new double[len + 1];
    beta[0] = 0.0;
    for (int i = 1; i <= len; ++i) beta[i] = 0.01;
  }
  else {
    beta[0] = 0.0;
    for (int i = 1; i <= len; ++i) fin>> beta[i];
  }
}

void DMSTransModel::write(std::ofstream fout) {
  fout<< len;

  fout.precision(10);
  fout.setf(0, std::ios::floatfield);

  if (beta == NULL) {
    for (int i = 1; i <= len; ++i) fout<< '\t'<< gamma[i];
  }
  else {
    for (int i = 1; i <= len; ++i) fout<< '\t'<< beta[i];
  }

  fout<< std::endl;
}
