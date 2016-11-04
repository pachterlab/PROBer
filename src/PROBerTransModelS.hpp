/* Copyright (c) 2016
   Bo Li (University of California, Berkeley)
   bli25@berkeley.edu

   This program is free software; you can redistribute it and/or
   modify it under the terms of the GNU General Public License as
   published by the Free Software Foundation; either version 3 of the
   License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
   General Public License for more details.   

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
   USA
*/

#ifndef PROBERTRANSMODELS_H_
#define PROBERTRANSMODELS_H_

#include <cmath>
#include <cassert>
#include <algorithm>
#include <string>
#include <fstream>
#include <vector>

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
class PROBerTransModelS {
public:
	/*
		@param   tid     transcript id (internal use)
		@param   name    transcript name
		@param   transcript_length  the length of this transcript
	 */
	PROBerTransModelS(int tid, const std::string& name = "", int transcript_length = -1);

	~PROBerTransModelS();

	/*
		@param   primer_length   the length of random primer
		@param   min_frag_len    minimum fragment length
		@param   max_frag_len    maximum fragment length
		@param   init_state      the initial state ((-) or (+) channel, if learning jointly or not)
		@comment: This function sets parameters shared by all transcripts for both simulation and learning, should be called before any PROBerTransModel object is created.
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
	static void setLearningRelatedParams(double gamma_init, double beta_init, double base, int read_length, bool isMAP, bool turnOnHidden);

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

	// In single transcript Model, state is not static anymore
	/*
		@return   current state
	*/
	int getState() { return state; }

	/*
		@comment: change channel
	 */
	void flipState() { state = state ^ 1; }

	/*
		@return   which channel we are dealing with (0, -; 1, +)
	 */
	int getChannel() const { return state & 1; }

	/*
		@return   if joint learning, true; otherwise, false
	*/
	bool isJoint() const { return state >= 2; }

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
		@return   N_obs, number of observed reads
	 */
	double getNobs() const { return N_obs[getChannel()]; }

	/*
		@param   channel   which channel to return
		@return the probability of passing the size selection step
	 */
	double getProbPass(int channel) const { return prob_pass[channel]; }

	/*                                                                                                                                                                                                         @param   channel   which channel to return                                                                                                                                                               @return  the log prob part of this transcript, 0 if this transcript is excluded                                                                                                                       */
	double getLogProbT(int channel) const { return log_prob[channel]; }

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

	// clear the alignments
	void clear() {
		N_obs[0] = N_obs[1] = 0.0;
		prob_pass[0] = prob_pass[1] = 1.0;
		log_prob[0] = log_prob[1] = 0.0;
		for (int i = 0; i < 2; ++i) {
			memset(starts[i], 0, sizeof(double) * (len + 1));
			memset(ends[i], 0, sizeof(double) * (len + 1));
			memset(ends_se[i], 0, sizeof(double) * (len + 1));
		}
	}

	/*
		@param   alignment   an in memory alignment belong to this transcript
		@return   true if the alignment is added, false otherwise
		@comment:  if the alignment's fragment length is not in [min_frag_len, max_frag_len] range, reject it
	 */
	bool addAlignment(InMemAlign* alignment) {
		int frag_len = alignment->fragment_length;
		int channel = getChannel();

		if (frag_len == 0) { // SE reads
			if (alignment->pos + min_alloc_len > len) return false;
			hasSE = true; // we have at least one SE read
			ends[channel][alignment->pos] += alignment->frac;
			ends_se[channel][alignment->pos] += alignment->frac;
			N_se[channel] += alignment->frac;
			N_obs[channel] += alignment->frac;
		}
		else { // PE reads
			frag_len -= primer_length;
			if (frag_len < min_frag_len || frag_len > max_frag_len || alignment->pos + frag_len > len)  return false;
			ends[channel][alignment->pos] += alignment->frac;
			starts[channel][alignment->pos + frag_len] += alignment->frac;
			N_obs[channel] += alignment->frac;
		}

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

	/*
		@param   choice    0, gamma; 1, beta
		@return  a copy of gamma or beta starting from position 1
	*/
	double* getCopy(int choice) {
		double* vec = new double[len];
		memcpy(vec, (choice == 0 ? gamma : beta) + 1, sizeof(double) * len);
		return vec;
	}

private:
	static const double INF; // Define exp(1000) as infinite to avoid the partial sum be -inf

	static int primer_length; // primer_length, the length of primers
	static int min_frag_len, max_frag_len; // min_frag_len and max_frag_len, the min and max fragment length (primer length excluded)

	static int init_state; // init state 

	static double gamma_init, beta_init;

	static int min_alloc_len; // minimum fragment length (primer length excluded) for allocating single end reads
	static bool isMAP; // if we should use MAP estimates instead of ML estimates
	static double base, dgamma, cgamma, dbeta, cbeta; // if MAP, gamma ~ Beta(dgamma + 1, cgamma + 1), beta ~ Beta(dbeta + 1, cbeta + 1); base = dgamma + cgamma = dbeta + cbeta

	static bool learning; // true if learning parameters, false if simulation

	static bool turnOnHidden; // if have hidden reads calculated

	int tid; // transcript id
	std::string name; // transcript name

	int state; // 0, learn/simulate gamma; 1, learn/simulate beta; 2, joint learning, gamma; 3, joint learning, beta

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

	/*
		comment: Auxiliary arrays below
	 */
	double *logsum; // logsum[i] = \sigma_{j=1}^{i} log(1-gamma[j]) if beta == NULL or \sigma_{j=1}^{i} log(1-gamma[j])(1-beta[j]). Thus a product from a to b is exp(logsum[b]-logsum[a-1]). 
	double *margin_prob; // margin_prob[i] = \sigma_{j = i + min_frag_len} ^ {i + max_frag_len} \prod_{k=i + min_frag_len + 1} ^{j} (1 - gamma[k]) * (beta == NULL ? 1.0 : (1 - beta[k]))

	int efflen2; // number of positions can generate a full length SE read
	double *margin_prob2; // not NULL only if min_alloc_len > min_frag_len, margin_prob2[i] = \sigma_{j = i + min_alloc_len} ^ {i + max_frag_len} \prod_{k = i + min_alloc_len + 1} ^ {j} (1 - gamma[k]) * (beta == NULL ? 1.0 : (1 - beta[k]))

	double *start2, *end2; // including hidden data, can be shared by a whole thread of transcripts

	double N_se[2]; // number of reads that we do not know their start positions
	double *starts[2], *ends[2], *ends_se[2];

	double log_prob[2]; // log prior + log lik of a transcript, lgamma function parts are excldued 

	double *cdf_end; // cumulative probabilities of having a read end at a particular position, only used for simulation


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
