#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <cassert>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <pthread.h>

#include "utils.h"
#include "my_assert.h"

#include "sampling.hpp"
#include "InMemoryStructs.hpp"
#include "PROBerTransModelS.hpp"

using namespace std;

bool verbose = true;

const int MAX_ROUND = 1000; // default maximum iterations
const double deltaChange = 5e-6; // default log probability change per read


struct Params {
	int no;
	Sampler *sampler;
	ifstream fin;
	vector<double*> estimates;
};

char *input_dir, *tname;
int target_id, target_length, num_trials, nThreads;

int primer_length, min_frag_len, max_frag_len, read_length;
double gamma_init, beta_init;

bool hasControl;

Params *paramsArray;
pthread_t *threads;
pthread_attr_t attr;
int rc;

bool hasSeed;
seedType seed;

EngineFactory factory;

void find_target(char* refName, string transcript_name) {
	char inpF[STRLEN];
	ifstream fin;
	string line;
	istringstream strin;

	string name;
	int id, len;

	target_id = target_length = 0;

	sprintf(inpF, "%s.translist", refName);
	fin.open(inpF);
	id = 0;
	while (getline(fin, line)) {
		strin.clear(); strin.str(line);
		strin>> name>> len;
		++id;

		if (name == transcript_name) {
			target_id = id;
			target_length = len;
			break;
		}
	}
	fin.close();

	assert(target_id > 0);
}

// assign threads
void init() {
	char inpF[STRLEN];

	paramsArray = new Params[nThreads];
	threads = new pthread_t[nThreads];

	hasSeed ? factory.init(seed) : factory.init();
	for (int i = 0; i < nThreads; ++i) {
		paramsArray[i].no = i;
		paramsArray[i].sampler = factory.new_sampler();		
		sprintf(inpF, "%s/gibbs_out_%d.txt", input_dir, i);
		paramsArray[i].fin.open(inpF);
		paramsArray[i].estimates.clear();
	}

	/* set thread attribute to be joinable */
	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

	if (verbose) { printf("Initialization finished!\n"); }
}

void runEM(PROBerTransModelS *model, int Nobs) {
	int channel;
	double prev_logprob, curr_logprob;
	double change;

	prev_logprob = curr_logprob = -1e300;

	// initialize
	model->init();
	
	// joint mode, EM
	channel = model->getChannel();
	assert(channel == 0);
	model->calcAuxiliaryArrays(channel);


	for (int ROUND = 1; ROUND <= MAX_ROUND; ++ROUND) {
		curr_logprob = model->getLogProbT(channel);
		model->EM_step(0.0);
		if (hasControl) {
			model->flipState();
			channel = model->getChannel();
			curr_logprob += model->getLogProbT(channel);
			model->EM_step(0.0);
			model->flipState();
		}
		channel = model->getChannel();

		change = (curr_logprob - prev_logprob) / Nobs;
		if (change <= deltaChange) break;

		prev_logprob = curr_logprob;
	}
}

void* bootstrap(void* arg) {
	Params *params = (Params*)arg;
	PROBerTransModelS *model;
	int tsize, csize;
	vector<InMemAlign> treatment, control;

	string line;
	int tid, pos, fragment_length;

	model = new PROBerTransModelS(target_id, tname, target_length);
	while (getline(params->fin, line)) {
		istringstream strin(line);
		while (strin>> tid>> pos>> fragment_length) {
			if (tid < 0) {
				tid = -tid;
				if (tid == target_id) control.push_back(InMemAlign(tid, pos, fragment_length, 0.0, 1.0));
			}
			else {
				if (tid == target_id) treatment.push_back(InMemAlign(tid, pos, fragment_length, 0.0, 1.0));
			}
		}

		tsize = treatment.size();
		csize = control.size();
		printf("tsize = %d, csize = %d\n", tsize, csize);
		
		for (int i = 0; i < num_trials; ++i) {
			model->clear();
			if (hasControl) {
				for (int j = 0; j < csize; ++j) 
					model->addAlignment(&control[int(params->sampler->random() * csize)]);
				model->flipState();
			}
			for (int j = 0; j < tsize; ++j)
				model->addAlignment(&treatment[int(params->sampler->random() * tsize)]);
			if (hasControl) model->flipState();

			runEM(model, (hasControl ? tsize + csize : tsize));
			params->estimates.push_back(model->getCopy(hasControl ? 1 : 0));
			printf("no = %d, i = %d\n", params->no, i);
		}
	}
	delete model;

	return NULL;
}

void output() {
	char outF[STRLEN];
	FILE *fo;
	int len = target_length - PROBerTransModelS::get_primer_length();
	
	int pos, vlen = 0;
	vector<double> values;

	for (int i = 0; i < nThreads; ++i)
		vlen += paramsArray[i].estimates.size();
	values.assign(vlen, 0.0);

	sprintf(outF, "%s/%s.txt", input_dir, tname);
	fo = fopen(outF, "w");
	for (int i = 0; i < len; ++i) {
		pos = 0;
		for (int j = 0; j < nThreads; ++j) 
			for (size_t k = 0; k < paramsArray[j].estimates.size(); ++k)
				values[pos++] = paramsArray[j].estimates[k][i];
		sort(values.begin(), values.end());
		for (int j = 0; j < vlen; ++j)	
			fprintf(fo, "%.6g%c", values[j], (j < vlen - 1 ? ' ' : '\n'));
	}
	fclose(fo);
}

void release() {
	pthread_attr_destroy(&attr);
	delete[] threads;

	for (int i = 0; i < nThreads; ++i) {
		paramsArray[i].fin.close();
		delete paramsArray[i].sampler;
		for (size_t k = 0; k < paramsArray[i].estimates.size(); ++k) 
			delete[] paramsArray[i].estimates[k];
	}
	delete[] paramsArray;			
}



int main(int argc, char* argv[]) {
	if (argc < 5) {
		printf("Usage: PROBer-bootstrap ref_name input_dir transcript_name num_trials" 
			" [--primer-length primer_length(default: 6)] [--size-selection-min min_frag_len(required)]"
			" [--size-selection-max max_frag_len(required)] [--read-length read_length]"
			" [--gamma-init gamma_init(default: 0.0001)] [--beta-init beta_init(default: 0.0001)]"
			" [-p number_of_threads] [--no-control] [--seed seed] [-q]\n");
		exit(-1);
	}

	input_dir = argv[2];
	tname = argv[3];
	num_trials = atoi(argv[4]);

	primer_length = 6;
	read_length = 0;
	min_frag_len = max_frag_len = -1;
	gamma_init = beta_init = 0.0001;

	nThreads = 1;
	hasControl = true;
	hasSeed = false;
	for (int i = 5; i < argc; ++i) {
		if (!strcmp(argv[i], "--primer-length")) primer_length = atoi(argv[i + 1]);
		if (!strcmp(argv[i], "--size-selection-min")) min_frag_len = atoi(argv[i + 1]);
		if (!strcmp(argv[i], "--size-selection-max")) max_frag_len = atoi(argv[i + 1]);
		if (!strcmp(argv[i], "--read-length")) read_length = atoi(argv[i + 1]);
		if (!strcmp(argv[i], "--gamma-init")) gamma_init = atof(argv[i + 1]);
		if (!strcmp(argv[i], "--beta-init")) beta_init = atof(argv[i + 1]);

		if (!strcmp(argv[i], "-p")) nThreads = atoi(argv[i + 1]);
		if (!strcmp(argv[i], "--no-control")) hasControl = false;
		if (!strcmp(argv[i], "--seed")) {
		  hasSeed = true;
		  int len = strlen(argv[i + 1]);
		  seed = 0;
		  for (int k = 0; k < len; ++k) seed = seed * 10 + (argv[i + 1][k] - '0');
		}
		if (!strcmp(argv[i], "-q")) verbose = false;
	} 

	find_target(argv[1], tname);
	assert(min_frag_len >= 0 && max_frag_len >= 0);
	PROBerTransModelS::setGlobalParams(primer_length, min_frag_len, max_frag_len, (hasControl ? 2 : 0));
	PROBerTransModelS::setLearningRelatedParams(gamma_init, beta_init, 1.0, read_length, true, false);
	init();
	for (int i = 0; i < nThreads; ++i) {
		rc = pthread_create(&threads[i], &attr, bootstrap, (void*)(&paramsArray[i]));
		pthread_assert(rc, "pthread_create", "Cannot create thread " + itos(i) + " (numbered from 0)!");
	}
	for (int i = 0; i < nThreads; ++i) {
		rc = pthread_join(threads[i], NULL);
		pthread_assert(rc, "pthread_join", "Cannot join thread " + itos(i) + " (numbered from 0)!");
	}
	printf("Bootstrapping is done!\n");
	output();
	release();

	return 0;
}
