#ifndef DMSTRANSMODEL_H_
#define DMSTRANSMODEL_H_

#include<cmath>
#include<cassert>
#include<algorithm>
#include<string>
#include<fstream>
#include<vector>

#include "utils.h"
#include "sampling.hpp"
#include "InMemoryStructs.hpp"

/*
  The coordinate system used outside is 0-based, starting from 5' end.
  The coordinate system used internally is 1-based, starting from 5' end:
     5'               3'
     0 - 1 2 ... n - n+1
   n = transcript length - primer length
 */
class DMSTransModel {
public:
  /*
    @param   learning    if we need to learn parameters
    @param   name    transcript name
    @param   transcript_length  the length of this transcript
   */
  DMSTransModel(bool learning, const std::string& name = "", int transcript_length = -1);

  ~DMSTransModel();

  /*
    @comment: This function sets parameters shared by all transcripts, should be called before any DMSTransModel object is created.
   */
  static void setGlobalParams(int primer_length, int min_frag_len, int max_frag_len, double gamma_init, double beta_init, int read_length);

  /*
    @return   primer length
   */
  static int get_primer_length() { return primer_length; }

  /*
    @return   minimum fragment length
   */
  static int get_minimum_fragment_length() { return min_frag_len + primer_length; }

  /*
    @return   maximum fragment length
   */
  static int get_maximum_fragment_length() { return max_frag_len + primer_length; }

  /*
    @return   transcript name
   */    
  const std::string& getName() const { return name; }

  /*
    @return   number of positions can estimate paramter: transcript length - primer length
   */
  int getLen() const { return len; }

  /*
    @return   number of positions can produce reads that pass the size selection step
   */
  int getEffLen() const { return efflen; }

  /*
    @return   N_obs, number of observed reads
   */
  double getNobs() const { return N_obs; }

  /*
    @comment: if this transcript has 0 reads, set its gamma/beta values to default
   */
  void setDefault() {
    if (!isZero(N_obs)) return;
    N_obs = 0.0;
    if (beta == NULL) memset(gamma, 0, sizeof(double) * (len + 1));
    else memset(beta, 0, sizeof(double) * (len + 1));
  }

  /*
    @param   pos     leftmost position in 5' end, 0-based  
    @return  the probability of generating a SE read end at pos
   */
  double getProb(int pos) const {
    int start_pos = pos + min_alloc_len;
    if (start_pos > len || pos < 0) return 0.0;
    double res = delta * (min_alloc_len == min_frag_len ? margin_prob[pos] : margin_prob2[pos]) * exp(logsum[start_pos] - logsum[pos]);
    if (pos > 0) res *= (beta == NULL ? gamma[pos] : (gamma[pos] + beta[pos] - gamma[pos] * beta[pos]));

    return res;
  }
  
  /*
    @param   pos      same as the above function
    @param   fragment_length     fragment length of the PE read
    @return  the probability of generating a PE read pair end at pos and has fragment length fragment_length
   */
  double getProb(int pos, int fragment_length) const {
    fragment_length -= primer_length;
    if (fragment_length < min_frag_len || fragment_length > max_frag_len) return 0.0;
    int start_pos = pos + fragment_length;
    if (start_pos > len || pos < 0) return 0.0;
    
    double res = delta * exp(logsum[start_pos] - logsum[pos]);
    if (pos > 0) res *= (beta == NULL ? gamma[pos] : (gamma[pos] + beta[pos] - gamma[pos] * beta[pos]));

    return res;
  }

  /*
    @param   alignment   an in memory alignment belong to this transcript
    @return   true if the alignment is added, false otherwise
    @comment:  if the alignment's fragment length is not in [min_frag_len, max_frag_len] range, reject it
   */
  bool addAlignment(InMemAlign* alignment) {
    int frag_len = alignment->fragment_length;
    if (frag_len == 0) { // SE reads
      if (alignment->pos + min_alloc_len > len) return false;
    }
    else { // PE reads
      frag_len -= primer_length;
      if (frag_len < min_frag_len || frag_len > max_frag_len || alignment->pos + frag_len > len)  return false;
    }

    alignments.push_back(alignment);

    return true;
  }

  /*
    @comment: This function calculate logsum and margin_prob and prob_pass, which are used to speed up the calculation
    @comment: It should be called before EM or getProb or get ProbPass etc. 
   */
  void calcAuxiliaryArrays();

  // Update counts information at each position
  void update();

  /*
    @return the probability of passing the size selection step
   */
  double getProbPass() const { return prob_pass; }

  /*
    @param   start2   auxiliary array used in EM
    @param   end2     auxiliary array used in EM
    @comment: start2 and end2 can be shared among multiple transcripts in a same thread. 
              This function tells this transcript where the start2 and end2 it should use.
   */
  void setStart2andEnd2(double* start2, double* end2) {
    this->start2 = start2;
    this->end2 = end2;
  }

  /*
    @param   N_tot   expected total counts for this transcript
    @comment: Run one iteration of EM algorithm  a single transcript
   */
  void EM(double N_tot);

  /*
    @param   fin   input stream
    @format:  name len [beta/gamma] ... 
   */
  void read(std::ifstream& fin);

  /*
    @param   fout output stream
    @format: the same as read
   */
  void write(std::ofstream& fout);

  /*
    @param   fc   output stream for c, the marking rate
    @param   fout output stream for freqs
    @format:  name c(rate of being marked) len freqs
   */
  void writeFreq(std::ofstream& fc, std::ofstream& fout);

  /*
    @comment: allocate memory for cdf_end, call calcAuxiliaryArrays() and calculate cdf_end values
   */
  void startSimulation();

  /*
    @param   sampler  a sampler used for sampling
    @param   pos      sampled 5' position (0-based)
    @param   fragment_length   sampled fragment length, with primer length considered
    @comment: If the probability of one fragment length is less than the MT19937 can generate, this fragment will never be generated. 
              If we use the naive simulation procedure, we may alleviate this problem.
   */
  void simulate(Sampler* sampler, int& pos, int& fragment_length);

  /*
    @comment: free memory allocated to cdf_end
   */
  void finishSimulation();

  double* getGamma() { return gamma; }
  double* getBeta() { return beta; }

private:
  static const double INF; // Define exp(1000) as infinite to avoid the partial sum be -inf

  static int primer_length; // primer_length, the length of primers
  static int min_frag_len, max_frag_len; // min_frag_len and max_frag_len, the min and max fragment length (primer length excluded)

  static double gamma_init, beta_init;

  static int min_alloc_len; // minimum fragment length (primer length excluded) for allocating single end reads


  std::string name; // transcript name
  bool learning; // if learn parameters
  int len; // len, number of position can learn parameters, transcript_length - primer_length
  int efflen; // efflen, number of positions can generate a valid fragment, len - min_frag_len + 1
  double N_obs; // Total number of observed counts
  double delta; // probability of priming from a particular position, delta = 1.0 / (len + 1)
  double prob_pass; // probability of generating a read that passes the size selection step
  double *gamma, *beta; // gamma, the vector of probability of drop-off at i (1-based); beta, the vector of probability of demtheylation at position i (1-based); 
  double *start, *end; // start, number of reads with first base after primer starting at a position; end, number of reads whose TF drops off at a position

  /*
    comment: Auxiliary arrays below
   */
  double *logsum; // logsum[i] = \sigma_{j=1}^{i} log(1-gamma[j]) if beta == NULL or \sigma_{j=1}^{i} log(1-gamma[j])(1-beta[j]). Thus a product from a to b is exp(logsum[b]-logsum[a-1]). 
  double *margin_prob; // margin_prob[i] = \sigma_{j = i + min_frag_len} ^ {i + max_frag_len} \prod_{k=i + min_frag_len + 1} ^{j} (1 - gamma[k]) * (beta == NULL ? 1.0 : (1 - beta[k]))

  int efflen2; // number of positions can generate a full length read
  double *margin_prob2; // not NULL only if min_alloc_len > min_frag_len, margin_prob2[i] = \sigma_{j = i + min_alloc_len} ^ {i + max_frag_len} \prod_{k = i + min_alloc_len + 1} ^ {j} (1 - gamma[k]) * (beta == NULL ? 1.0 : (1 - beta[k]))

  double *start2, *end2; // including hidden data, can be shared by a whole thread of transcripts

  double *cdf_end; // cumulative probabilities of having a read end at a particular position, only used for simulation

  std::vector<InMemAlign*> alignments; // In memory alignments used for update
};

#endif
