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

#ifndef FRAGFragLenDist_H_
#define FRAGFragLenDist_H_

#include <cassert>
#include <fstream>
#include <vector>

class FragLenDist {
public:
  FragLenDist(int maxL = 1000);

  int getMinL() const { return lb; }
  int getMaxL() const { return ub; }

  double getProb(int len) const { return pmf[len - lb]; }

  // For estimating FragLenDist, distribute multi-mapping reads evenly
  // I have looked at one real eCLIP data set, the distribution of unique reads only and all reads (with multi reads evenly distributed) are very similar
  void update(int len, double frac = 1.0) {
    if (len > ub) { ub = len; pmf.resize(len, 0.0); }
    pmf[len - lb] += frac;
  }

  void finish();

  void read(std::ifstream& fin);
  void write(std::ofstream& fout);

private:
  int lb, ub, span; // [lb, ub], span = ub - lb + 1
  std::vector<double> pmf; // probability mass function, cumulative density function, and counts of noise reads
};

#endif 
