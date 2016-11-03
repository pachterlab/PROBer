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

#include "utils.h"
#include "NoiseProfile.hpp"

NoiseProfile::NoiseProfile(bool hasCount) {
	memset(p, 0, sizeof(p));
	c = pc = NULL;
	if (hasCount) { 
		c = new double[NCODES];
		memset(c, 0, sizeof(double) * NCODES);
	}
}

NoiseProfile::~NoiseProfile() {
	if (c != NULL) delete[] c;
}

void NoiseProfile::writeC(std::ofstream& fout) {
	fout<< NCODES<< std::endl;
	fout.precision(0);
	fout.setf(std::ios::fixed, std::ios::floatfield);
	for (int i = 0; i < NCODES; ++i) fout<< c[i]<< '\t';
	fout<< std::endl;
}

// read sufficient statistics from unalignable reads
void NoiseProfile::readC(std::ifstream& fin) {
	int tmp_ncodes;
	assert((fin>> tmp_ncodes) && (tmp_ncodes == NCODES));
	for (int i = 0; i < NCODES; ++i) assert(fin>> c[i]);
}

void NoiseProfile::calcInitParams() {
	double sum = 0.0;

	// one pseudo count for each base 
	for (int i = 0; i < NCODES; ++i) sum += (1.0 + c[i]);
	for (int i = 0; i < NCODES; ++i) p[i] = (1.0 + c[i]) / sum;
}

void NoiseProfile::init() {
	memset(p, 0, sizeof(p));
}

void NoiseProfile::collect(const NoiseProfile* o) {
	for (int i = 0; i < NCODES; ++i)
		p[i] += o->p[i];
}

void NoiseProfile::finish() {
	double sum = 0.0;
	
	for (int i = 0; i < NCODES; ++i) sum += (p[i] + c[i]);
	if (isZero(sum)) memset(p, 0, sizeof(p)); 
	else for (int i = 0; i < NCODES; ++i) p[i] = (p[i] + c[i]) / sum;
}

double NoiseProfile::calcLogP() {
	double logp = 0.0;
	for (int i = 0; i < NCODES; ++i) 
		if (c[i] > 0.0) logp += c[i] * log(p[i]);
	return logp;
}

void NoiseProfile::read(std::ifstream& fin) {
	std::string line;
	while (getline(fin, line)) {
		if (line.substr(0, 13) == "#NoiseProfile") break;
	}
	assert(fin.good());

	int tmp_ncodes;
	assert((fin>> tmp_ncodes) && (tmp_ncodes == NCODES));

	for (int i = 0; i < NCODES; ++i)
		assert(fin>> p[i]);

	getline(fin, line);
}

void NoiseProfile::write(std::ofstream& fout) {
	fout<< "#NoiseProfile, format: NCODES; P_noise for each base"<< std::endl;
	fout<< NCODES<< std::endl;

	for (int i = 0; i < NCODES - 1; ++i) fout<< p[i]<< '\t';
	fout<< p[NCODES - 1]<< std::endl<< std::endl;
}

void NoiseProfile::startSimulation() {
	pc = new double[NCODES];
	
	memcpy(pc, p, sizeof(p));
	for (int i = 1; i < NCODES; ++i) pc[i] += pc[i - 1];
}

void NoiseProfile::finishSimulation() {
	delete[] pc;
}
