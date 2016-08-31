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

#include <cstring>
#include <cassert>
#include <string>
#include <fstream>

#include "utils.h"
#include "Profile.hpp"

Profile::Profile(int maxL) {
  proLen = maxL;
  size = proLen * NCODES * NCODES;
  p = new double[proLen][NCODES][NCODES];
  memset(p, 0, sizeof(double) * size);
  
  //set initial parameters
  int N = NCODES - 1;
  double probN = 1e-5, portionC = 0.99; //portionC, among ACGT, the portion of probability mass the correct base takes
  double probC, probO;
  
  for (int i = 0; i < proLen; ++i) {
    for (int j = 0; j < NCODES - 1; ++j) {
      p[i][j][N] = probN;
      probC = portionC * (1.0 - probN);
      probO = (1.0 - portionC) / (NCODES - 2) * (1.0 - probN);
      
      for (int k = 0; k < NCODES - 1; ++k) {
	p[i][j][k] = (j == k ? probC : probO);
      }
    }
    p[i][N][N] = probN;
    for (int k = 0; k < NCODES - 1; ++k)
      p[i][N][k] = (1.0 - probN) / (NCODES - 1);
  }

  pc = NULL;
}

Profile::~Profile() { 
  delete[] p;
}

void Profile::init() {
  memset(p, 0, sizeof(double) * size);
}

void Profile::collect(const Profile* o) {
  for (int i = 0; i < proLen; ++i)
    for (int j = 0; j < NCODES; ++j)
      for (int k = 0; k < NCODES; ++k)
	p[i][j][k] += o->p[i][j][k];
}

void Profile::finish() {
  double sum;
  
  for (int i = 0; i < proLen; ++i) {
    for (int j = 0; j < NCODES; ++j) {
      sum = 0.0;
      for (int k = 0; k < NCODES; ++k) sum += p[i][j][k];
      if (isZero(sum)) memset(p[i][j], 0, sizeof(double) * NCODES);
      else for (int k = 0; k < NCODES; ++k) p[i][j][k] /= sum;
    }
  }
}

void Profile::read(std::ifstream& fin) {
  std::string line;
  while (getline(fin, line)) {
    if (line.substr(0, 8) == "#Profile") break;
  }
  assert(fin.good());

  int tmp_prolen, tmp_ncodes;
  fin>> tmp_prolen>> tmp_ncodes;
  assert(fin.good() && (tmp_ncodes == NCODES));

  if (tmp_prolen != proLen) {
    delete[] p;
    proLen = tmp_prolen;
    size = proLen * NCODES * NCODES;
    p = new double[proLen][NCODES][NCODES];
    memset(p, 0, sizeof(double) * size);
  }
  
  for (int i = 0; i < proLen; ++i)
    for (int j = 0; j < NCODES; ++j)
      for (int k = 0; k < NCODES; ++k)
	assert(fin>> p[i][j][k]);

  getline(fin, line);
}

void Profile::write(std::ofstream& fout) {
  fout<< "#Profile, format: proLen NCODES; P[POS][REF_BASE][OBSERVED_BASE], proLen blocks separated by a blank line, each block contains NCODES lines"<< std::endl;
  fout<< proLen<< '\t'<< NCODES<< std::endl;

  for (int i = 0; i < proLen; ++i) {
    for (int j = 0; j < NCODES; ++j) {
      for (int k = 0; k < NCODES - 1; ++k)
	fout<< p[i][j][k]<< '\t';
      fout<< p[i][j][NCODES - 1]<< std::endl;
    }
    fout<< std::endl;
  }
}

void Profile::startSimulation() {
  pc = new double[proLen][NCODES][NCODES];

  for (int i = 0; i < proLen; ++i) {
    for (int j = 0; j < NCODES; ++j)
      for (int k = 0; k < NCODES; ++k) {
	pc[i][j][k] = p[i][j][k];
	if (k > 0) pc[i][j][k] += pc[i][j][k - 1];
      }

    // In case one (Pos, Ref_base) combination is not seen, sharing information from other combinations  
    double cp_sum, cp_d, cp_n;
    double p_d, p_o, p_n;

    cp_sum = cp_d = cp_n = 0.0;
    for (int j = 0; j < NCODES - 1; ++j) {
      cp_sum += pc[i][j][NCODES - 1];
      cp_d += p[i][j][j];
      cp_n += p[i][j][NCODES - 1];
    }
    
    if (isZero(cp_sum)) {
      p_n = 1e-5;
      p_d = 0.99 * (1.0 - p_n);
      p_o = (1.0 - p_d - p_n) / (NCODES - 2);
    }
    else {
      p_d = cp_d / cp_sum;
      p_n = cp_n / cp_sum;
      p_o = (1.0 - p_d - p_n) / (NCODES - 2);
    }

    // Check if (Pos, j) has no probability for j != N
    for (int j = 0; j < NCODES - 1; ++j) {
      if (!isZero(pc[i][j][NCODES - 1])) continue;
      
      for (int k = 0; k < NCODES; ++k) {
	if (k == j) pc[i][j][k] = p_d;
	else if (k == NCODES - 1) pc[i][j][k] = p_n;
	else pc[i][j][k] = p_o;
	if (k > 0) pc[i][j][k] += pc[i][j][k - 1];
      }
    }
    
    // Check if (Pos, N) has no probability
    if (isZero(pc[i][NCODES - 1][NCODES - 1])) {
      p_o = (1.0 - p_n) / (NCODES - 1);
      for (int k = 0; k < NCODES; ++k) {
	pc[i][NCODES - 1][k] = (k < NCODES - 1 ? p_o : p_n);
	if (k > 0) pc[i][NCODES - 1][k] += pc[i][NCODES - 1][k - 1];
      }
    } 
  }  
}

void Profile::finishSimulation() {
  delete[] pc;
}

