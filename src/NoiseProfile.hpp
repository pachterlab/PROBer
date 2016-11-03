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

#ifndef NOISEPROFILE_H_
#define NOISEPROFILE_H_

#include <cstring>
#include <string>
#include <fstream>

#include "utils.h"
#include "sampling.hpp"
#include "SEQstring.hpp"

class NoiseProfile {
public:
	NoiseProfile(bool hasCount = false);
	~NoiseProfile();
	
	void updateC(const SEQstring& seq) {
		int len = seq.getLen();
		for (int i = 0; i < len; ++i) {
			++c[seq.baseCodeAt(i)];
		}
	}

	void writeC(std::ofstream& fout);
	void readC(std::ifstream& fin);

	void calcInitParams();

	double getProb(const SEQstring& seq) {
		double prob = 1.0;
		int len = seq.getLen();
		
		for (int i = 0; i < len; ++i) {
			prob *= p[seq.baseCodeAt(i)];
		}
		
		return prob;
	}
	
	void update(const SEQstring& seq, double frac) {
		int len = seq.getLen();
		for (int i = 0; i < len; ++i) {
			p[seq.baseCodeAt(i)] += frac;
		}
	}

	void init();
	void collect(const NoiseProfile* o);
	void finish();
	
	double calcLogP();

	void read(std::ifstream& fin);
	void write(std::ofstream& fout);

	void simulate(Sampler* sampler, int len, std::string& readseq) {
		readseq.assign(len, 0);
		for (int i = 0; i < len; ++i) {
			readseq[i] = code2base[sampler->sample(pc, NCODES)];
		}
	}
	
	void startSimulation();
	void finishSimulation();
	
private:
	static const int NCODES = 5;

	double p[NCODES];
	double *c; // counts in N0

	double *pc; // for simulation
};

#endif /* NOISEPROFILE_H_ */
