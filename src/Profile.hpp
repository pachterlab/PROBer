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

#ifndef PROFILE_H_
#define PROFILE_H_

#include <fstream>

#include "utils.h"
#include "sampling.hpp"

class Profile {
public:
	Profile(int maxL = 1000);
	~Profile();

	double getProb(int pos, int ref_base, int read_base) {
		return p[pos][ref_base][read_base];
	}

	void update(int pos, int ref_base, int read_base, double frac) {
		p[pos][ref_base][read_base] += frac;
	}

	void init();
	void collect(const Profile* o);
	void finish(); // No pseudo-count here. However, if we assume Illumina platform and the read length difference is due to trimming, it should be fine.

	void read(std::ifstream& fin);
	void write(std::ofstream& fout);

	char simulate(Sampler* sampler, int pos, int ref_base) {
		return code2base[sampler->sample(pc[pos][ref_base], NCODES)];
	}

	void startSimulation();
	void finishSimulation();
	
private:
	static const int NCODES = 5;
	
	int proLen; // profile length
	int size; // # of items in p;
	double (*p)[NCODES][NCODES]; //profile matrices
	
	double (*pc)[NCODES][NCODES]; // for simulation
};

#endif /* PROFILE_H_ */
