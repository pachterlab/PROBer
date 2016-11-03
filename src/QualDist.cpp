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
#include "QualDist.hpp"

QualDist::QualDist() {
	memset(p_init, 0, sizeof(p_init));
	memset(p_tran, 0, sizeof(p_tran));
}

double QualDist::finish() {
	double logp;
	double sum, value;
	
	logp = 0.0;

	sum = 0.0;
	for (int i = 0; i < SIZE; ++i) sum += p_init[i];
	assert(!isZero(sum));
	for (int i = 0; i < SIZE; ++i) {
		value = p_init[i] / sum;
		if (p_init[i] > 0.0) logp += p_init[i] * log(value);
		p_init[i] = value;
	}
	
	for (int i = 0; i < SIZE; ++i) {
		sum = 0.0;
		for (int j = 0; j < SIZE; ++j) sum += p_tran[i][j];
		if (isZero(sum)) memset(p_tran[i], 0, sizeof(double) * SIZE);
		else for (int j = 0; j < SIZE; ++j) {
	value = p_tran[i][j] / sum;
	if (p_tran[i][j] > 0.0) logp += p_tran[i][j] * log(value);
	p_tran[i][j] = value;
			}
	}

	return logp;
}

void QualDist::read(std::ifstream& fin) {
	std::string line;
	while (getline(fin, line)) {
		if (line.substr(0, 9) == "#QualDist") break;
	}
	assert(fin.good());

	int tmp_size;
	assert((fin>> tmp_size) && (tmp_size == SIZE));

	for (int i = 0; i < SIZE; ++i) assert(fin>> p_init[i]);
	for (int i = 0; i < SIZE; ++i) 
		for (int j = 0; j < SIZE; ++j) assert(fin>> p_tran[i][j]);

	getline(fin, line);
}

void QualDist::write(std::ofstream& fout) {
	fout<< "#QualDist, format: SIZE; P_init; P_tran, SIZE lines"<< std::endl;
	fout<< SIZE<< std::endl;
	
	for (int i = 0; i < SIZE - 1; ++i) fout<< p_init[i]<< '\t';
	fout<< p_init[SIZE - 1]<< std::endl;
	for (int i = 0; i < SIZE; ++i) {
		for (int j = 0; j < SIZE -1 ; ++j) fout<< p_tran[i][j]<< '\t';
		fout<< p_tran[i][SIZE - 1]<< std::endl;
	}
	fout<< std::endl;
}

void QualDist::startSimulation() {
	qc_init = new double[SIZE];
	qc_trans = new double[SIZE][SIZE];
	
	memcpy(qc_init, p_init, sizeof(p_init));
	for (int i = 1; i < SIZE; ++i) qc_init[i] += qc_init[i - 1];

	memcpy(qc_trans, p_tran, sizeof(p_tran));
	for (int i = 0; i < SIZE; ++i)
		for (int j = 1; j < SIZE; ++j) 
			qc_trans[i][j] += qc_trans[i][j - 1];
}

void QualDist::finishSimulation() {
	delete[] qc_init;
	delete[] qc_trans;
}
