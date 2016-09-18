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

#include <ctime>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <cassert>
#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <unordered_map>

#include "htslib/sam.h"

#include "utils.h"
#include "my_assert.h"
#include "sampling.hpp"

#include "BamAlignment.hpp"
#include "AlignmentGroup.hpp"
#include "SamParser.hpp"
#include "BamWriter.hpp"

#include "PROBerReadModel_iCLIP.hpp"

using namespace std;

bool verbose = true; // define verbose

const int MAXALIGN = 105;

// position, key for map
struct KeyType {
	int cid, dir, pos; // dir, 0, +, 1, -

	KeyType() : cid(0), dir(0), pos(0) {}
	
	KeyType(int cid, char dir, int pos) : cid(cid), dir(dir == '+' ? 0 : 1), pos(pos) {}
	
	bool operator== (const KeyType& o) const {
		return (cid == o.cid) && (dir == o.dir) && (pos == o.pos);
	}
};

struct MyHash {
	size_t hash_combine(size_t seed, size_t value) const {
		return seed ^ (value + 0x9e3779b9 + (seed<< 6) + (seed>> 2));
	}

	size_t operator() (const KeyType& key) const {
		return hash_combine(hash_combine(key.cid, key.dir), key.pos);
	}
};

typedef unordered_map<string, int> DictType;

typedef unordered_map<KeyType, double, MyHash> HashType;
typedef HashType::iterator HashIterType;



char *sampleName, *statName;
char *alignF;

PROBerReadModel_iCLIP* model = NULL;

const bam_hdr_t* header = NULL;
BamAlignment* ba;
AlignmentGroup ag;
SamParser* parser = NULL;
BamWriter* writer = NULL;

DictType chr2cid;
HashType weightTable;
HashIterType it;

int mate; // The mate carry binding information, default is 1, which assumes iCLIP protocol; 2 if eCLIP
bool isNaive;

seedType seed;
Sampler* sampler = NULL;

void init() {
	char inpF[STRLEN];

	// Build chr2cid table
	parser = new SamParser(alignF);
	header = parser->getHeader();
	for (int i = 0; i < header->n_targets; ++i)
		chr2cid[string(header->target_name[i])] = i;

	// Load sequencing error model
	if (!isNaive) {
		sprintf(inpF, "%s.model", statName);
		model = new PROBerReadModel_iCLIP();
		model->read(inpF);
	}

	// Load site info
	string chr;
	char dir;
	int pos, nuniq;
	double nmulti, tot;

	int cnt = 0;

	sprintf(inpF, "%s.site_info", sampleName);
	ifstream fin(inpF);
	while (fin>> chr>> dir>> pos>> nuniq>> nmulti) {
		tot = nuniq + nmulti;
		if (tot > 0.0) weightTable[KeyType(chr2cid[chr], dir, pos)] = tot;

		++cnt;
		if (verbose && (cnt % 5000000 == 0)) cout<< "Loaded "<< cnt<< " sites."<< endl; 
	}
	fin.close();
	if (verbose) cout<< "init is finished."<< endl;
}

void sampleMultiReads() {
	int s;
	char outF[STRLEN];
	vector<double> fracs(MAXALIGN, 0.0);
	double conprbs[MAXALIGN];

	sampler = new Sampler(seed);
	sprintf(outF, "%s.sampled.bam", sampleName);
	writer = new BamWriter(outF, header, "Sampled iCLIP/eCLIP BAM");

	READ_INT_TYPE cnt = 0, N11 = 0, N12 = 0;

	do {
		ag.clear();
		
		while (parser->next(ag)) {
			if (!ag.isFiltered() && ag.isAligned()) {
				s = ag.size();
				if (s == 1) { ++N11; writer->write(ag); }
				else {
					fracs.assign(s, 1.0);
					if (!isNaive) {
						model->calcProbs(ag, conprbs);
						for (int i = 0; i < s; ++i) fracs[i] *= conprbs[i];
					}
					for (int i = 0; i < s; ++i) {
						ba = ag.getAlignment(i);
						it = weightTable.find(KeyType(ba->getTid(), ba->getMateDir(mate), ba->getCrosslinkSite(mate)));
						fracs[i] *= (it == weightTable.end() ? 0.0 : it->second);
						if (i > 0) fracs[i] += fracs[i - 1];
					}
					if (fracs[s - 1] > 0.0) {
						ba = ag.getAlignment(sampler->sample(fracs, s));
						writer->write(*ba);
						++N12;
					}
				}
			}

			++cnt;
			if (verbose && (cnt % 1000000 == 0)) cout<< cnt<< " reads are processed!"<< endl;
		}
		
		if (verbose) cout<< alignF<< " is processed."<< endl;

		delete parser;
		parser = NULL;

		alignF = strtok(NULL, ",");
		if (alignF != NULL) parser = new SamParser(alignF);
	} while (parser != NULL);

	if (verbose) cout<< "N11 = "<< N11<< ", N12 = "<< N12<< "."<< endl;

	delete writer;
	delete sampler;
}

int main(int argc, char* argv[]) {
	if (argc < 4) { 
		printf("PROBer-sample-iCLIP sampleName statName alignFileList [--eCLIP] [--naive] [--seed seed] [-q]\n");
		exit(-1);
	}

	sampleName = argv[1];
	statName = argv[2];
	alignF = strtok(argv[3], ",");


	mate = 1;
	isNaive = false;
	seed = time(NULL);
	for (int i = 4; i < argc; ++i) {
		if (!strcmp(argv[i], "-q")) verbose = false;
		if (!strcmp(argv[i], "--eCLIP")) mate = 2;
		if (!strcmp(argv[i], "--naive")) isNaive = true;
		if (!strcmp(argv[i], "--seed")) {
			seed = 0;
			int len = strlen(argv[i + 1]);
			for (int j = 0; j < len; ++j) seed = seed * 10 + (argv[i + 1][j] - '0');
		}
	}

	init();
	sampleMultiReads();

	return 0;
}
