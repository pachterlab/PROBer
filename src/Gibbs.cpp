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
#include <cstring>
#include <cstdlib>
#include <cassert>
#include <fstream>
#include <sstream>
#include <vector>
#include <pthread.h>

#include "utils.h"
#include "my_assert.h"

#include "sampling.hpp"
#include "InMemoryStructs.hpp"

using namespace std;

bool verbose = true;

struct Params {
	int no, nsamples;
	FILE *fo;
	Sampler *sampler;
};

int nThreads;
char *outdir, *controlF;

int M;
READ_INT_TYPE N0, N1;

int BURNIN, NSAMPLES, GAP;

double pseudoC;
vector<HIT_INT_TYPE> s;
vector<InMemAlign> hits;

bool quiet;

Params *paramsArray;
pthread_t *threads;
pthread_attr_t attr;
int rc;

bool hasSeed;
seedType seed;

EngineFactory factory;

void load_data(char* inpF, bool isControl) {
	ifstream fin;
	string line;
	int tmp_N1, tmp_N0, tmp_M;

	fin.open(inpF);
	general_assert(fin.is_open(), "Cannot open " + cstrtos(inpF) + "!");
	fin>> tmp_N1>> tmp_N0>> tmp_M;
	if (!isControl) N1 = tmp_N1, N0 = tmp_N0, M = tmp_M;
	else N1 += tmp_N1, N0 += tmp_N0, assert(M == tmp_M); 
	getline(fin, line);

	while (getline(fin, line)) {
		istringstream strin(line);

		int size;
		int tid, pos, fragment_length;
		double conprb;

		strin>> size>> conprb;
		hits.push_back(InMemAlign(0, 0, 0, conprb));
		for (int i = 0; i < size; ++i) {
			strin>> tid>> pos>> fragment_length>> conprb;
			hits.push_back(InMemAlign(isControl ? -tid : tid, pos, fragment_length, conprb));
		}
		s.push_back(hits.size());
	}

	fin.close();

	if (verbose) { printf("Loading %s is finished.\n", inpF); }
}


// assign threads
void init() {
	int quotient, left;
	char outF[STRLEN];

	quotient = NSAMPLES / nThreads;
	left = NSAMPLES % nThreads;

	paramsArray = new Params[nThreads];
	threads = new pthread_t[nThreads];

	hasSeed ? factory.init(seed) : factory.init();
	for (int i = 0; i < nThreads; ++i) {
		paramsArray[i].no = i;

		paramsArray[i].nsamples = quotient;
		if (i < left) ++paramsArray[i].nsamples;

		sprintf(outF, "%s/gibbs_out_%d.txt", outdir, i);
		paramsArray[i].fo = fopen(outF, "w");

		paramsArray[i].sampler = factory.new_sampler();		
	}

	/* set thread attribute to be joinable */
	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

	if (verbose) { printf("Initialization finished!\n"); }
}

void writeGibbsOut(FILE* fo, vector<int>& z) {
	int hid;
	for (int i = 0; i < N1; ++i) {
		hid = s[i] + z[i];
		if (hits[hid].tid != 0) 
			fprintf(fo, "%d %d %d ", hits[hid].tid, hits[hid].pos, hits[hid].fragment_length);
	}
	fprintf(fo, "\n");
}

void* Gibbs(void* arg) {
	int CHAINLEN;
	HIT_INT_TYPE len, fr, to;
	Params *params = (Params*)arg;

	vector<int> z(N1, 0), counts(M + 1, 0);
	vector<double> arr;

	// generate initial state
	counts[0] = N0;
	for (READ_INT_TYPE i = 0; i < N1; ++i) {
		fr = s[i]; to = s[i + 1];
		len = to - fr;
		arr.assign(len, 0);
		for (HIT_INT_TYPE j = fr; j < to; ++j) {
			arr[j - fr] = hits[j].conprb;
			if (j > fr) arr[j - fr] += arr[j - fr - 1];  // cumulative
		}
		z[i] = abs(hits[fr + params->sampler->sample(arr, len)].tid);
		++counts[z[i]];
	}

	// Gibbs sampling
	CHAINLEN = 1 + (params->nsamples - 1) * GAP;
	for (int ROUND = 1; ROUND <= BURNIN + CHAINLEN; ++ROUND) {

		for (READ_INT_TYPE i = 0; i < N1; ++i) {
			--counts[z[i]];
			fr = s[i]; to = s[i + 1]; len = to - fr;
			arr.assign(len, 0);
			for (HIT_INT_TYPE j = fr; j < to; j++) {
				arr[j - fr] = (counts[abs(hits[j].tid)] + pseudoC) * hits[j].conprb;
				if (j > fr) arr[j - fr] += arr[j - fr - 1]; //cumulative
			}
			z[i] = abs(hits[fr + params->sampler->sample(arr, len)].tid);
			++counts[z[i]];
		}

		if ((ROUND > BURNIN) && ((ROUND - BURNIN - 1) % GAP == 0))
			writeGibbsOut(params->fo, z);

		if (verbose && ROUND % 100 == 0) { printf("Thread %d, ROUND %d is finished!\n", params->no, ROUND); }
	}

	return NULL;
}

void release() {
	/* destroy attribute */
	pthread_attr_destroy(&attr);
	delete[] threads;

	for (int i = 0; i < nThreads; ++i) {
		fclose(paramsArray[i].fo);
		delete paramsArray[i].sampler;
	}
	delete[] paramsArray;
}

int main(int argc, char* argv[]) {
	if (argc < 6) {
		printf("Usage: PROBer-run-gibbs plus_channel_input output_dir BURNIN NSAMPLE GAP [--minus minus_channel_input] [-p #Threads] [--seed seed] [--pseudo-count pseudo_count] [-q]\n");
		exit(-1);
	}

	outdir = argv[2];

	BURNIN = atoi(argv[3]);
	NSAMPLES = atoi(argv[4]);
	GAP = atoi(argv[5]);

	controlF = NULL;
	nThreads = 1;
	hasSeed = false;
	pseudoC = 1.0;
	quiet = false;

	for (int i = 6; i < argc; ++i) {
		if (!strcmp(argv[i], "--minus")) controlF = argv[i + 1];
		if (!strcmp(argv[i], "-p")) nThreads = atoi(argv[i + 1]);
		if (!strcmp(argv[i], "--seed")) {
		  hasSeed = true;
		  int len = strlen(argv[i + 1]);
		  seed = 0;
		  for (int k = 0; k < len; ++k) seed = seed * 10 + (argv[i + 1][k] - '0');
		}
		if (!strcmp(argv[i], "--pseudo-count")) pseudoC = atof(argv[i + 1]);
		if (!strcmp(argv[i], "-q")) quiet = true;
	}
	verbose = !quiet;

	assert(NSAMPLES > 1); // Otherwise, we cannot calculate posterior variance

	if (nThreads > NSAMPLES) {
		nThreads = NSAMPLES;
		fprintf(stderr, "Warning: Number of samples is less than number of threads! Change the number of threads to %d!\n", nThreads);
	}

	s.clear(); s.push_back(0); hits.clear();
	load_data(argv[1], false);
	if (controlF != NULL) load_data(controlF, true);
	assert(N1 == s.size() - 1);

	if (verbose) printf("Gibbs started!\n");
	system(("mkdir -p " + string(outdir)).c_str());

	init();
	for (int i = 0; i < nThreads; ++i) {
		rc = pthread_create(&threads[i], &attr, Gibbs, (void*)(&paramsArray[i]));
		pthread_assert(rc, "pthread_create", "Cannot create thread " + itos(i) + " (numbered from 0)!");
	}
	for (int i = 0; i < nThreads; ++i) {
		rc = pthread_join(threads[i], NULL);
		pthread_assert(rc, "pthread_join", "Cannot join thread " + itos(i) + " (numbered from 0)!");
	}
	release();

	if (verbose) printf("Gibbs finished!\n");

	return 0;
}
