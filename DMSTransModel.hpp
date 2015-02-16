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
    @param   tid     transcript id (internal use)
    @param   name    transcript name
    @param   transcript_length  the length of this transcript
   */
  DMSTransModel(int tid, const std::string& name = "", int transcript_length = -1);

  ~DMSTransModel();

  /*
    @param   primer_length   the length of random primer
    @param   min_frag_len    minimum fragment length
    @param   max_frag_len    maximum fragment length
    @param   init_state      the initial state ((-) or (+) channel, if learning jointly or not)
    @comment: This function sets parameters shared by all transcripts for both simulation and learning, should be called before any DMSTransModel object is created.
   */
  static void setGlobalParams(int primer_length, int min_frag_len, int max_frag_len, int init_state);
  
  /*
    @param   gamma_init   initial value for gamma
    @param   beta_init    initial value for beta
    @param   base         base = alpha + beta - 2 for Beta(alpha, beta) distribution
    @param   read_length  the length of a single-end read
    @param   isMAP        if we want MAP estimates
    @comment: This function sets learning related parameters shared by all transcripts, calling this function means that we want to learn parameters from data.
  */
  static void setLearningRelatedParams(double gamma_init, double beta_init, double base, int read_length, bool isMAP);

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
    @return   if learning or simulation
   */
  static bool isLearning() { return learning; }

  /*
    @return   current state
  */
  static int getState() { return state; }

  /*
    @comment: change channel
   */
  static void flipState() { state = state ^ 1; }

  /*
    @return   which channel we are dealing with (0, -; 1, +)
   */
  static int getChannel() { return state & 1; }

  /*
    @return   if joint learning, true; otherwise, false
  */
  static bool isJoint() { return state >= 2; }

  /*
    @return   transcript id
  */
  int getTid() const { return tid; }

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
    @param   channel   which channel to look at
    @return   number of alignments this transcript has for channel
  */
  HIT_INT_TYPE getNumAlignments(int channel) const {
    return alignmentsArr[channel].size();
  }

  /*
    @comment: check if this transcript can be excluded from learning procedure due to either no available position or no alignments. Must call after alignments are processed!
  */
  bool isExcluded() const {
    return efflen <= 0 || (!isJoint() && alignmentsArr[getChannel()].size() == 0) || (isJoint() && alignmentsArr[0].size() == 0 && alignmentsArr[1].size() == 0);
  }

  /*
    @return   N_obs, number of observed reads
   */
  double getNobs() const { return N_obs[getChannel()]; }

  /*
    @return the probability of passing the size selection step
   */
  double getProbPass() const { return prob_pass[getChannel()]; }

  /*
    @param   pos     leftmost position in 5' end, 0-based  
    @return  the probability of generating a SE read end at pos
   */
  double getProb(int pos) const {
    int start_pos = pos + min_alloc_len;
    if (start_pos > len || pos < 0) return 0.0;
    double res = delta * (min_alloc_len == min_frag_len ? margin_prob[pos] : margin_prob2[pos]) * exp(logsum[start_pos] - logsum[pos]);
    if (pos > 0) res *= (getChannel() == 0 ? gamma[pos] : (gamma[pos] + beta[pos] - gamma[pos] * beta[pos]));

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
    if (pos > 0) res *= (getChannel() == 0 ? gamma[pos] : (gamma[pos] + beta[pos] - gamma[pos] * beta[pos]));

    return res;
  }

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
    @comment: if this transcript has 0 reads, set its gamma/beta values to default
   */

  /*
    @param   alignment   an in memory alignment belong to this transcript
    @return   true if the alignment is added, false otherwise
    @comment:  if the alignment's fragment length is not in [min_frag_len, max_frag_len] range, reject it
   */
  bool addAlignment(InMemAlign* alignment) {
    int frag_len = alignment->fragment_length;
    if (frag_len == 0) { // SE reads
      if (alignment->pos + min_alloc_len > len) return false;
      hasSE = true; // we have at least one SE read
    }
    else { // PE reads
      frag_len -= primer_length;
      if (frag_len < min_frag_len || frag_len > max_frag_len || alignment->pos + frag_len > len)  return false;
    }

    alignmentsArr[getChannel()].push_back(alignment);

    return true;
  }

  /*
    @comment: Initialize related data members to prepare this transcript for parameter esitmation. Call only after all alignments are added.
   */
  void init();

  /*
    @param   channel   which channel we should calculate for
    @comment: This function calculate logsum and margin_prob and prob_pass, which are used to speed up the calculation
    @comment: It should be called before EM or getProb or get ProbPass etc. 
   */
  void calcAuxiliaryArrays(int channel);

  // Update counts information at each position
  void update();

  /*
    @param   N_tot   expected total counts for this transcript
    @comment: Run one iteration of EM algorithm for a single transcript
   */
  void EM_step(double N_tot);

  /*
    @param   fin   input stream
    @param   channel   which channel
    @format:  name len [beta/gamma] ... 
   */
  void read(std::ifstream& fin, int channel);

  /*
    @param   fout     output stream
    @param   channel  which channel 
    @format: the same as read
   */
  void write(std::ofstream& fout, int channel);

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
  static int state; // 0, learn/simulate gamma; 1, learn/simulate beta; 2, joint learning, gamma; 3, joint learning, beta

  static double gamma_init, beta_init;

  static int min_alloc_len; // minimum fragment length (primer length excluded) for allocating single end reads
  static bool isMAP; // if we should use MAP estimates instead of ML estimates
  static double base, dgamma, cgamma, dbeta, cbeta; // if MAP, gamma ~ Beta(dgamma + 1, cgamma + 1), beta ~ Beta(dbeta + 1, cbeta + 1); base = dgamma + cgamma = dbeta + cbeta

  static bool learning; // true if learning parameters, false if simulation

  int tid; // transcript id
  std::string name; // transcript name

  int len; // len, number of position can learn parameters, transcript_length - primer_length
  int efflen; // efflen, number of positions can generate a valid fragment, len - min_frag_len + 1
  double delta; // probability of priming from a particular position, delta = 1.0 / (len + 1)
  double N_obs[2]; // Total number of observed counts
  double prob_pass[2]; // probability of generating a read that passes the size selection step
  double *gamma, *beta; // gamma, the vector of probability of drop-off at i (1-based); beta, the vector of probability of demtheylation at position i (1-based); 
  double *start, *end; // start, number of reads with first base after primer starting at a position; end, number of reads whose TF drops off at a position
  double *dcm, *ccm; // drop-off counts and covering counts for (-) channel
  double *end_se; // number of SE reads end at a position

  bool hasSE; // if this transcript has SE reads, which means we do not know its start, from either of the channels
  double N_se; // number of reads that we do not know their start positions

  /*
    comment: Auxiliary arrays below
   */
  double *logsum; // logsum[i] = \sigma_{j=1}^{i} log(1-gamma[j]) if beta == NULL or \sigma_{j=1}^{i} log(1-gamma[j])(1-beta[j]). Thus a product from a to b is exp(logsum[b]-logsum[a-1]). 
  double *margin_prob; // margin_prob[i] = \sigma_{j = i + min_frag_len} ^ {i + max_frag_len} \prod_{k=i + min_frag_len + 1} ^{j} (1 - gamma[k]) * (beta == NULL ? 1.0 : (1 - beta[k]))

  int efflen2; // number of positions can generate a full length SE read
  double *margin_prob2; // not NULL only if min_alloc_len > min_frag_len, margin_prob2[i] = \sigma_{j = i + min_alloc_len} ^ {i + max_frag_len} \prod_{k = i + min_alloc_len + 1} ^ {j} (1 - gamma[k]) * (beta == NULL ? 1.0 : (1 - beta[k]))

  double *start2, *end2; // including hidden data, can be shared by a whole thread of transcripts

  double *cdf_end; // cumulative probabilities of having a read end at a particular position, only used for simulation

  std::vector<InMemAlign*> alignmentsArr[2]; // In memory alignments used for update from (-) and (+) channels

  /*
    @param   beta   beta value at a position, this is the to-be-estimated parameter
    @param   gamma  gamma value at the same position, which is assumed known
    @param   dc     drop-off counts at the position
    @param   cc     covering counts at the position
    @comment: This function estimate MAP beta by solving an quadratic equation. It is only used when we estimate gamma and beta separately.
   */
  void solveQuadratic1(double& beta, double gamma, double dc, double cc);

  /*
    @param   gamma   gamma to be estimated
    @param   beta    beta to be estimated
    @param   dcm     drop-off counts at (-) channel
    @param   ccm     covering counts at (-) channel
    @param   dcp     drop-off counts at (+) channel
    @param   ccp     covering counts at (+) channel
   */
  void solveQuadratic2(double& gamma, double& beta, double dcm, double ccm, double dcp, double ccp);
};

#endif
