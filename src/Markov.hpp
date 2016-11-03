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

#ifndef MARKOV_H_
#define MARKOV_H_

#include <vector>
#include <fstream>

#include "utils.h"
#include "sampling.hpp"

class Markov {
public:
	Markov();

	double getProb(char a) { return P_start[chr2state[a]]; }
	double getProb(char a, char b) { return P_trans[chr2state[a]][chr2state[b]]; }
	
	double getIBaseProb(int code) { return probI[code]; }

	void update(char a, double frac) { P_start[chr2state[a]] += frac; }
	void update(char a, char b, double frac) { P_trans[chr2state[a]][chr2state[b]] += frac; }

	void updateIBase(int code, double frac) { probI[code] += frac; }

	void init();
	void collect(const Markov* o);
	void finish();

	void read(std::ifstream& fin);
	void write(std::ofstream& fout);

	char simulate(Sampler *sampler) {
		switch(sampler->sample(P_start_sim, NSTATES)) {
		case M : return 'M';
		case I : return 'I';
		case D : return 'D';
		default : assert(false);
		}
	}

	char simulate(Sampler *sampler, char a) {
		switch(sampler->sample(P_trans_sim[chr2state[a]], NSTATES)) {
		case M : return 'M';
		case I : return 'I';
		case D : return 'D';
		default : assert(false);
		}
	}

	char simulateIBase(Sampler *sampler) {
		return code2base[sampler->sample(probI_sim, NCODES)];
	}

	void startSimulation();
	void finishSimulation();

private:
	static const int NSTATES = 3;
	static const int M = 0;
	static const int I = 1;
	static const int D = 2;

	static const std::vector<int> chr2state;
	static std::vector<int> init_chr2state();

	double P_start[NSTATES]; // model the probability of a read starts from a state
	double P_trans[NSTATES][NSTATES];

	static const int NCODES = 5; // ACGTN
	
	double probI[NCODES]; // The probability of generating a base given the state is I

	// for simulation
	double *P_start_sim;
	double (*P_trans_sim)[NSTATES];
	double *probI_sim;
};

#endif 
