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
#include <cassert>
#include <string>
#include <fstream>

#include "MateLenDist.hpp"

MateLenDist::MateLenDist() {
	lb = ub = span = 0;
	logp = 0.0;
	
	pmf.clear();
	cdf.clear();
	noise_counts.clear();
}

void MateLenDist::finish() {
	double sum = 0.0;

	assert(ub > 0);
	for (lb = 1; pmf[lb] == 0.0; ++lb);
	for (int i = lb; i <= ub; ++i) {
		sum += pmf[i];
		pmf[i - lb] = pmf[i];
	}
	span = ub - lb + 1;
	pmf.resize(span);
	for (int i = 0; i < span; ++i) pmf[i] /= sum;

	cdf.assign(span, 0.0);
	for (int i = 0; i < span; ++i) cdf[i] = pmf[i] + (i > 0 ? cdf[i - 1] : 0.0);

	// Calculate logp
	logp = 0.0;
	for (int i = 0; i < span; ++i) 
		if (noise_counts[i + lb] > 0.0) logp += noise_counts[i + lb] * log(pmf[i]);
}

void MateLenDist::read(std::ifstream& fin) {
	std::string line;
	while (getline(fin, line)) {
		if (line.substr(0, 12) == "#MateLenDist") break;
	}
	assert(fin.good());

	assert(fin>> lb>> ub>> span);
	pmf.resize(span, 0.0);
	cdf.resize(span, 0.0);
	for (int i = 0; i < span; ++i) {
		assert(fin>> pmf[i]);
		cdf[i] = pmf[i];
		if (i > 0) cdf[i] += cdf[i - 1];
	}

	getline(fin, line);
}

void MateLenDist::write(std::ofstream& fout) {
	fout<< "#MateLenDist, format: lb ub span; [lb, ub], span = ub - lb + 1, probability mass function values"<< std::endl;
	fout<< lb<< '\t'<< ub<< '\t'<< span<< std::endl;  
	for (int i = 0; i < span - 1; ++i) fout<< pmf[i]<< '\t';
	fout<< pmf[span - 1]<< std::endl<< std::endl;
}
