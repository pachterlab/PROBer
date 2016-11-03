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

#include <cmath>
#include <cstring>
#include <cassert>
#include <string>
#include <fstream>
#include <algorithm>

#include "utils.h"
#include "sampling.hpp"
#include "PROBerTransModelS.hpp"

const double PROBerTransModelS::INF = 1000.0;

int PROBerTransModelS::primer_length = 6; // default, 6bp
int PROBerTransModelS::min_frag_len;
int PROBerTransModelS::max_frag_len;

double PROBerTransModelS::gamma_init;
double PROBerTransModelS::beta_init;

double PROBerTransModelS::base;
double PROBerTransModelS::cgamma;
double PROBerTransModelS::dgamma;
double PROBerTransModelS::cbeta;
double PROBerTransModelS::dbeta;

int PROBerTransModelS::min_alloc_len;
bool PROBerTransModelS::isMAP = true; // default is true

bool PROBerTransModelS::learning = false; // default is simulation
int PROBerTransModelS::init_state = 0;

bool PROBerTransModelS::turnOnHidden = false; // default is not to turn on


void PROBerTransModelS::setGlobalParams(int primer_length, int min_frag_len, int max_frag_len, int init_state) { 
	assert(primer_length <= min_frag_len && min_frag_len <= max_frag_len);
	PROBerTransModelS::primer_length = primer_length;
	PROBerTransModelS::min_frag_len = min_frag_len - primer_length;
	PROBerTransModelS::max_frag_len = max_frag_len - primer_length;
	PROBerTransModelS::init_state = init_state;
}

void PROBerTransModelS::setLearningRelatedParams(double gamma_init, double beta_init, double base, int read_length, bool isMAP, bool turnOnHidden) {
	learning = true;

	PROBerTransModelS::gamma_init = gamma_init;
	PROBerTransModelS::beta_init = beta_init;
	PROBerTransModelS::base = base;
	PROBerTransModelS::min_alloc_len = (read_length < min_frag_len ? min_frag_len : read_length) - primer_length;
	PROBerTransModelS::isMAP = isMAP;
	PROBerTransModelS::turnOnHidden = turnOnHidden;

	if (isMAP) {
		dgamma = gamma_init * base;
		cgamma = base - dgamma;
		dbeta = beta_init * base;
		cbeta = base - dbeta;
	}
}


PROBerTransModelS::PROBerTransModelS(int tid, const std::string& name, int transcript_length) : tid(tid), name(name) {  
	state = PROBerTransModelS::init_state; // initialize the state of this transcript

	gamma = beta = NULL;
	start = end = NULL;
	dcm = ccm = NULL;
	end_se = NULL;

	logsum = margin_prob = NULL;
	margin_prob2 = NULL;

	start2 = end2 = NULL;

	len = efflen = -1; 
	efflen2 = -1;
	N_obs[0] = N_obs[1] = 0.0;
	prob_pass[0] = prob_pass[1] = 1.0; // In case no alignments, unobserved read counts is 0

	delta = 0.0;

	hasSE = false;

	N_se[0] = N_se[1] = 0.0;
	for (int i = 0; i < 2; ++i) {
		starts[i] = NULL;
		ends[i] = NULL;
		ends_se[i] = NULL;
	}

	log_prob[0] = log_prob[1] = 0.0;
	cdf_end = NULL;

	if (!learning) return;

	len = transcript_length - primer_length;
	efflen = len - min_frag_len + 1;
	delta = 1.0 / (len + (primer_length > 0 ? 1.0 : 0.0));  

	// If a transcript is excluded from analysis, all its gamma/beta values become 0
	gamma = new double[len + 1];
	memset(gamma, 0, sizeof(double) * (len + 1));

	if (getState() > 0) {
		beta = new double[len + 1];
		memset(beta, 0, sizeof(double) * (len + 1));
	}

	for (int i = 0; i < 2; ++i) {
		starts[i] = new double[len + 1];
		ends[i] = new double[len + 1];
		ends_se[i] = new double[len + 1];
		memset(starts[i], 0, sizeof(double) * (len + 1));
		memset(ends[i], 0, sizeof(double) * (len + 1));
		memset(ends_se[i], 0, sizeof(double) * (len + 1));
	}  
}

PROBerTransModelS::~PROBerTransModelS() {
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

	if (start2 != NULL) delete[] start2;
	if (end2 != NULL) delete[] end2;

	for (int i = 0; i < 2; ++i) {
		if (starts[i] != NULL) delete[] starts[i];
		if (ends[i] != NULL) delete[] ends[i];
		if (ends_se[i] != NULL) delete[] ends_se[i];
	}
}

void PROBerTransModelS::init() {
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

	// start2 and end2
	start2 = new double[len + 1];
	end2 = new double[len + 1];
	memset(start2, 0, sizeof(double) * (len + 1));
	memset(end2, 0, sizeof(double) * (len + 1));

	// Auxiliary arrays
	logsum = new double[len + 1];
	assert(efflen > 0);
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

void PROBerTransModelS::calcAuxiliaryArrays(int channel) {
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

	// Calculate log prob for this channel
	log_prob[channel] = 0.0;
	// prior
	if (isMAP) {
		if (channel == 0)
			for (int i = 1; i <= len; ++i) log_prob[channel] += dgamma * log(gamma[i]) + cgamma * log (1.0 - gamma[i]);
		else
			for (int i = 1; i <= len; ++i) log_prob[channel] += dbeta * log(beta[i]) + cbeta * log(1.0 - beta[i]);
	}
	// SE reads
	if (N_se[channel] > 0.0) 
		for (int i = 0; i < efflen2; ++i) 
			if (ends_se[channel][i] > 0.0) {
	value = delta * (min_alloc_len == min_frag_len ? margin_prob[i] : margin_prob2[i]) * exp(logsum[i + min_alloc_len] - logsum[i]);
				if (i > 0) value *= (channel == 0 ? gamma[i] : (gamma[i] + beta[i] - gamma[i] * beta[i]));
	assert(value > 0.0 && value < 1.0);
	log_prob[channel] += ends_se[channel][i] * log(value);
			}
	// PE reads
	if (N_obs[channel] > N_se[channel]) {
		log_prob[channel] += (N_obs[channel] - N_se[channel]) * log(delta); // log probability of starting paired-end reads or full fragments

		double ncover = std::max(ends[channel][0] - ends_se[channel][0] - starts[channel][0], 0.0);
		double ndrop;

		for (int i = 1; i <= len; ++i) {
			ndrop = std::max(ends[channel][i] - ends_se[channel][i], 0.0);
			if (ndrop > 0.0) log_prob[channel] += ndrop * log(channel == 0 ? gamma[i] : gamma[i] + beta[i] - gamma[i] * beta[i]);
			if (ncover > 0.0) log_prob[channel] += ncover * log(channel == 0 ? 1.0 - gamma[i] : (1.0 - gamma[i]) * (1.0 - beta[i]));
			ncover = std::max(ncover + ndrop - starts[channel][i], 0.0);
		}
	}
	// normalizing factor
	if (turnOnHidden && N_obs[channel] > 0.0) {
		assert(prob_pass[channel] > 0.0 && prob_pass[channel] < 1.0);
		log_prob[channel] -= N_obs[channel] * log(prob_pass[channel]);
	}
}

inline void PROBerTransModelS::solveQuadratic1(double& beta, double gamma, double dc, double cc) {
	double a = (1.0 - gamma) * (cbeta + cc + dbeta + dc);
	double b = ((cbeta + cc + 2.0 * dbeta + dc) * gamma - (dc + dbeta)) / a;
	double c = (-dbeta * gamma) / a;
	double sqt_delta = sqrt(b * b - 4.0 * c);

	//  if (sqt_delta <= fabs(b)) { printf("gamma = %.10g, dc = %.10g, cc = %.10g, a = %.10g, b = %.10g, c = %.10g, sqt_delta = %.10g\n", gamma, dc, cc, a, b, c, sqt_delta); }

	assert(sqt_delta > fabs(b));
	beta = (-b + sqt_delta) / 2.0;
	assert(beta > 0.0 && beta < 1.0);
}

inline void PROBerTransModelS::solveQuadratic2(double& gamma, double& beta, double dcm, double ccm, double dcp, double ccp) {
	double common_factor = cgamma + ccm - cbeta - dbeta;
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

void PROBerTransModelS::EM_step(double N_tot) {
	int channel = getChannel();

	int max_end_i;
	double psum, value;

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
		// copy vectors
		memcpy(start, starts[channel], sizeof(double) * (len + 1));
		memcpy(end, ends[channel], sizeof(double) * (len + 1));
		if (!isZero(N_se[channel])) memcpy(end_se, ends_se[channel], sizeof(double) * (len + 1));

		//E step, if we have reads that do not know their start positions, infer start from end
		if (!isZero(N_se[channel])) {
			double prev, curr;
			double *mp = NULL;
			int effl;

			if (min_alloc_len > min_frag_len) { mp = margin_prob2; effl = efflen2; }
			else { mp = margin_prob; effl = efflen; }
			assert(effl > 0);

			curr = (end_se[0] > 0.0 && mp[0] > 0.0) ? end_se[0] / mp[0] : 0.0;
			start[min_alloc_len] += curr;
			for (int i = 1, pos = min_alloc_len + 1; i < effl; ++i, ++pos) {
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

		if (turnOnHidden) {
			
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

		}
		else {
			memset(start2, 0, sizeof(double) * (len + 1));
			memset(end2, 0, sizeof(double) * (len + 1));
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
		assert(gamma[i] > 0.0 && gamma[i] < 1.0);
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

void PROBerTransModelS::read(std::ifstream& fin, int channel) {
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

void PROBerTransModelS::write(std::ofstream& fout, int channel) {
	fout<< name<< '\t'<< len;

	if (channel == 0) {
		for (int i = 1; i <= len; ++i) fout<< '\t'<< gamma[i];
	}
	else {
		for (int i = 1; i <= len; ++i) fout<< '\t'<< beta[i];
	}

	fout<< std::endl;
}

void PROBerTransModelS::writeFreq(std::ofstream& fc, std::ofstream& fout) {
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

void PROBerTransModelS::startSimulation() {
	if (efflen <= 0) return;
	cdf_end = new double[efflen];

	calcAuxiliaryArrays(getChannel());
 
	cdf_end[0] = exp(logsum[min_frag_len] - logsum[0]) * margin_prob[0];
	for (int i = 1; i < efflen; ++i) {
		cdf_end[i] = (getChannel() == 0 ? gamma[i] : (gamma[i] + beta[i] - gamma[i] * beta[i])) * exp(logsum[i + min_frag_len] - logsum[i]) * margin_prob[i];
		cdf_end[i] += cdf_end[i - 1];
	}
}

void PROBerTransModelS::simulate(Sampler* sampler, int& pos, int& fragment_length) {
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

void PROBerTransModelS::finishSimulation() {
	if (cdf_end != NULL) delete[] cdf_end;
	cdf_end = NULL;
}
