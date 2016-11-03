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

#include <cstdio>
#include <cassert>
#include <string>
#include <fstream>
#include <vector>

#include "my_assert.h"
#include "RefSeq.hpp"
#include "Refs.hpp"

Refs::Refs() {
	M = 0;
	seqs.assign(1, NULL);
}

Refs::~Refs() {
	for (int i = 1; i <= M; ++i) delete seqs[i];
}

void Refs::readFrom(char* inpF) {
	std::ifstream fin(inpF);
	RefSeq* seq;

	general_assert(fin.is_open(), "Cannot open " + cstrtos(inpF) + "! It may not exist.");
	
	M = 0;
	seqs.assign(1, NULL);
	seq = new RefSeq();
	
	while (seq->read(fin)) {
		seqs.push_back(seq);
		++M;
		seq = new RefSeq();
	}
	delete seq;
	
	fin.close();

	assert(M + 1 == (int)seqs.size());

	if (verbose) { printf("Refs.readFrom finished!\n"); }
}

void Refs::writeTo(char* outF) {
	std::ofstream fout(outF);
	for (int i = 1; i <= M; ++i) 
		seqs[i]->write(fout);
	fout.close();
	if (verbose) { printf("Refs.writeTo finished!\n"); }
}

void Refs::writeTransListTo(char* outF) {
	std::ofstream fout(outF);
	for (int i = 1; i <= M; ++i)
		fout<< seqs[i]->getName()<< '\t'<< seqs[i]->getLen()<< std::endl;
	fout.close();
	if (verbose) { printf("Refs.writeTransListTo finished!\n"); }
}
