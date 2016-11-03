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

#ifndef MATELENDIST_H_
#define MATELENDIST_H_

#include <cassert>
#include <fstream>
#include <vector>

#include "sampling.hpp"

class MateLenDist {
public:
	MateLenDist();

	int getMinL() const { return lb; }
	int getMaxL() const { return ub; }

	double getProb(int len, int upper_bound = MAXV) const {
		assert(len >= lb && len <= ub);
		return (len != upper_bound ? pmf[len - lb] : pmf[len - lb] + (cdf[span - 1] - cdf[len - lb]));
	}
		
	void update(int len, bool is_noise = false) {
		if (len > ub) { ub = len; pmf.resize(ub + 1, 0.0); noise_counts.resize(ub + 1, 0.0); }
		++pmf[len];
		if (is_noise) ++noise_counts[len];
	}

	void finish();

	// Get log probability for unalignable reads, call after finish()
	double getLogP() const { return logp; }

	void read(std::ifstream& fin);
	void write(std::ofstream& fout);

	int simulate(Sampler* sampler, int upper_bound = MAXV) {
		if (upper_bound < lb) return upper_bound;
		int len = lb + sampler->sample(cdf, span);
		if (upper_bound < len) len = upper_bound;
		return len;
	}

private:
	static const int MAXV = 2147483647; // the maximum int value, used as an upper bound on mate length

	int lb, ub, span; // [lb, ub], span = ub - lb + 1
	double logp; // probability of generating the mate lengths for unalignable reads
	std::vector<double> pmf, cdf, noise_counts; // probability mass function, cumulative density function, and counts of noise reads
};

#endif 
