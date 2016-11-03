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

#ifndef QPROFILE_H_
#define QPROFILE_H_

#include <cassert>
#include <fstream>

#include "utils.h"
#include "sampling.hpp"

class QProfile {
public:
	QProfile();
 
	// qual starts from 0, 33 is already deducted
	double getProb(int qual, int ref_base, int read_base) {
		return p[qual][ref_base][read_base];
	}

	void update(int qual, int ref_base, int read_base, double frac) {
		p[qual][ref_base][read_base] += frac;
	}

	void init();
	void collect(const QProfile* o);
	void finish();
		
	void read(std::ifstream& fin);
	void write(std::ofstream& fout);
 
	char simulate(Sampler* sampler, int qual, int ref_base) {
		return code2base[sampler->sample(pc[qual][ref_base], NCODES)];
	}
	
	void startSimulation();
	void finishSimulation();
	
private:
	static const int NCODES = 5; // number of possible codes
	static const int SIZE = 100;
	
	double p[SIZE][NCODES][NCODES]; // p[q][r][c] = p(c|r,q)
	
	double (*pc)[NCODES][NCODES]; // for simulation
};

#endif /* QPROFILE_H_ */
