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

#include <cmath>
#include <ctime>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <cassert>
#include <vector>
#include <string>
#include <iostream>
#include <pthread.h>

#include "htslib/sam.h"

#include "utils.h"
#include "my_assert.h"

#include "Refs.hpp"
#include "Transcripts.hpp"

#include "BamAlignment.hpp"
#include "AlignmentGroup.hpp"
#include "SamParser.hpp"
#include "BamWriter.hpp"

#include "PROBerWholeModel.hpp"
#include "PROBerReadModel.hpp"
#include "InMemoryStructs.hpp"

using namespace std;

bool verbose = true; // define verbose


const int MAX_ROUND = 1000; // default maximum iterations
const double deltaChange = 5e-6; // default log probability change per read

// Parameter struct to pass parameters to each subprocess
struct InMemParams {
	int no; // thread number

	PROBerWholeModel *whole_model;
	PROBerReadModel *read_model;

	PROBerReadModel *estimator; // slave model that is used to estimate model parameters
	InMemChunk *chunk; // A chunk of memory to record all in-memory information for reads and alignments associated with this thread

	double count0; // sum of noise read fractions
	double loglik; // log likelihood

	InMemParams(int no, PROBerWholeModel* whole_model, PROBerReadModel* read_model, READ_INT_TYPE nreads, HIT_INT_TYPE nlines) {
		this->no = no;
		this->whole_model = whole_model;
		this->read_model = read_model;
		estimator = NULL;
		chunk = new InMemChunk(nreads, nlines);
		count0 = loglik = 0.0;
	}

	~InMemParams() {
		delete chunk;
		if (estimator != NULL) delete estimator;
	}
};

bool keepGoing;

int M; // Number of transcripts
int N0[2], N_eff[2]; // Number of unalignable reads, number of effective reads (unaligned + aligned)
double count0[2], logprob[2]; // Used in EM algorithm, number of unalignable reads and log probability for each channel

int model_type; 
int num_threads;
int read_length;
bool isMAP;
bool has_control;

char refName[STRLEN], sampleName[STRLEN], imdName[STRLEN], statName[STRLEN], channel[STRLEN];

Refs refs;
Transcripts transcripts;

PROBerWholeModel *whole_model;
PROBerReadModel *read_models[2];

bool needCalcConPrb, updateReadModel;
vector<InMemParams*> paramsVecs[2];
vector<pthread_t> threads;
pthread_attr_t attr;
int rc;

bool output_bam, output_logMAP;

bam_hdr_t *hdr;

// Preprocess reads and alignments
void preprocessAlignments(int channel) {
	char bamF[STRLEN], partitionF[STRLEN];
	SamParser *parser = NULL;
	AlignmentGroup ag;

	if (verbose) { printf("Begin to preprocess BAM files for channel %s!\n", channelStr[channel]); }

	// Preprocess unalignable reads
	sprintf(bamF, "%s_%s_N0.bam", imdName, channelStr[channel]);
	parser = new SamParser(bamF);

	N0[channel] = 0;
	while (parser->next(ag)) {
		read_models[channel]->update_preprocess(ag, false);
		++N0[channel];
	}
	if (hdr == NULL) hdr = parser->pass_header();
	delete parser;
	if (verbose) { printf("Unalignable reads are preprocessed!\n"); }

	// Preprocess alignable reads
	sprintf(partitionF, "%s_%s.partition", imdName, channelStr[channel]);
	ifstream fin(partitionF);
	assert(fin.is_open());

	int id;
	READ_INT_TYPE rid, nreads;
	HIT_INT_TYPE nlines;

	InMemAlignG *a_read = NULL;
	InMemAlign *aligns = NULL;

	int cnt = 0;
	InMemChunk *chunk = NULL;

	bool is_paired;
	int seqlen;

	N_eff[channel] = N0[channel];
	paramsVecs[channel].assign(num_threads, NULL);
	for (int i = 0; i < num_threads; ++i) {
		fin>> id>> nreads>> nlines;
		assert(id == i);
		N_eff[channel] += nreads;
		paramsVecs[channel][i] = new InMemParams(i, whole_model, read_models[channel], nreads, nlines);

		sprintf(bamF, "%s_%s_%d.bam", imdName, channelStr[channel], i);
		parser = new SamParser(bamF, hdr);
		rid = 0;
		ag.clear();
		while (parser->next(ag)) {
			is_paired = ag.isPaired();
			seqlen = !is_paired ? ag.getSeqLength() : 0; 

			assert(paramsVecs[channel][i]->chunk->next(a_read, aligns));
			a_read->size = ag.size();
			
			for (int j = 0; j < ag.size(); ++j) {
				BamAlignment *ba = ag.getAlignment(j);
				aligns[j].tid = transcripts.getInternalSid(ba->getTid());
				aligns[j].pos = ba->getLeftMostPos();

				if (is_paired) aligns[j].fragment_length = ba->getInsertSize();
				else if (seqlen < read_length) aligns[j].fragment_length = ba->getAlignedLength();
				else aligns[j].fragment_length = 0;

				aligns[j].frac = 1.0 / ag.size();
			}
			whole_model->addAlignments(a_read, aligns);
			read_models[channel]->update_preprocess(ag, true);
			++rid;

			if (verbose && (rid % 1000000 == 0)) cout<< "Loaded "<< rid<< " reads!"<< endl;
		}
		assert(rid == nreads);

		delete parser;
		if (verbose) { printf("Thread %d's data is preprocessed!\n", i); }

		chunk = paramsVecs[channel][i]->chunk;
		for (HIT_INT_TYPE j = 0; j < nlines; ++j) 
			if (chunk->aligns[j].conprb == -1.0) ++cnt;
	}

	if (verbose) { printf("There are %d alignments filtered!\n", cnt); }

	for (int i = 0; i < num_threads; ++i) 
		paramsVecs[channel][i]->estimator = new PROBerReadModel(read_models[channel]);
	read_models[channel]->finish_preprocess();
	
	if (verbose) { printf("Bam preprocessing is done for channel %s!\n", channelStr[channel]); }

	// change to the other state for further analysis
	if (has_control) whole_model->flipState();
}

void init() {
	char refF[STRLEN], tiF[STRLEN];
	char configF[STRLEN];

	// Load references
	sprintf(refF, "%s.transcripts.fa", refName);
	refs.readFrom(refF);
	M = refs.getM();
	
	sprintf(tiF, "%s.ti", refName);
	transcripts.readFrom(tiF);
	char imd_name[STRLEN];
	sprintf(imd_name, "%s_plus", imdName);
	transcripts.buildMappings(imd_name);

	// Create PROBerWholeModel
	sprintf(configF, "%s.config", imdName);
	whole_model = new PROBerWholeModel(configF, (has_control ? 2 : 0), has_control, &transcripts, num_threads, read_length, isMAP);

	// Create PROBerReadModels
	read_models[0] = has_control ? new PROBerReadModel(model_type, &refs, read_length) : NULL;
	read_models[1] = new PROBerReadModel(model_type, &refs, read_length);

	memset(N0, 0, sizeof(N0));
	memset(N_eff, 0, sizeof(N_eff));

	memset(count0, 0, sizeof(count0));
	memset(logprob, 0, sizeof(logprob));
	
	hdr = NULL;
	// Preprocess data for (-)
	if (has_control) preprocessAlignments(0);
	// Preprocess data for (+)
	preprocessAlignments(1);

	threads.assign(num_threads, pthread_t());
	/* set thread attribute to be joinable */
	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

	if (verbose) { printf("Preprocess data is finished!\n"); }
}

void* E_STEP(void* arg) {
	InMemParams *params = (InMemParams*)arg;
	PROBerWholeModel *whole_model = params->whole_model;
	PROBerReadModel *read_model = params->read_model;
	PROBerReadModel *estimator = params->estimator;
	InMemChunk *chunk = params->chunk;

	SamParser *parser = NULL;
	AlignmentGroup ag;

	params->count0 = 0.0;
	params->loglik = 0.0;

	chunk->reset();

	if (needCalcConPrb || updateReadModel) {
		char bamF[STRLEN];
		sprintf(bamF, "%s_%s_%d.bam", imdName, whole_model->get_channel_string(whole_model->getChannel()), params->no);
		parser = new SamParser(bamF, hdr); 
	}
	if (updateReadModel) estimator->init();

	READ_INT_TYPE nreads = chunk->nreads;
	int size;
	double sum, noise_frac;
	InMemAlignG *a_read = NULL;
	InMemAlign *aligns = NULL;

	for (READ_INT_TYPE i = 0; i < nreads; ++i) {
		if (needCalcConPrb || updateReadModel) {
			assert(parser->next(ag));
		}

		assert(chunk->next(a_read, aligns));

		if (needCalcConPrb) read_model->setConProbs(a_read, aligns, ag);

		size = a_read->size;
		sum = noise_frac = whole_model->getProb(0) * a_read->noise_conprb;
		for (int j = 0; j < size; ++j) {
			if (aligns[j].conprb > 0.0) aligns[j].frac = whole_model->getProb(aligns[j].tid, aligns[j].pos, aligns[j].fragment_length) * aligns[j].conprb;
			else aligns[j].frac = 0.0;
			sum += aligns[j].frac;
		}
		assert(sum > 0.0);

		params->loglik += log(sum);
		noise_frac /= sum;
		params->count0 += noise_frac;
		for (int j = 0; j < size; ++j) aligns[j].frac /= sum;

		if (updateReadModel) estimator->update(a_read, aligns, ag, noise_frac);    
	}

	if (parser != NULL) delete parser;

	return NULL;
}

inline bool needUpdateReadModel(int ROUND) {
	return ROUND <= 10;
}

void one_EM_iteration(int channel, int ROUND) {
	// init
	if (ROUND == 1) whole_model->init();

	logprob[channel] = (isMAP ? whole_model->getLogPrior() : 0.0);

	// E step
	for (int i = 0; i < num_threads; ++i) {
		rc = pthread_create(&threads[i], &attr, E_STEP, (void*)paramsVecs[channel][i]);
		pthread_assert(rc, "pthread_create", "Cannot create thread " + itos(i) + " (numbered from 0) at EM ROUND " + itos(ROUND) + " for " + channelStr[channel] + " channel!");
	}
		
	for (int i = 0; i < num_threads; ++i) {
		rc = pthread_join(threads[i], NULL);
		pthread_assert(rc, "pthread_join", "Cannot join thread " + itos(i) + " (numbered from 0) at EM ROUND " + itos(ROUND) + " for " + channelStr[channel] + " channel!");
	}

	count0[channel] = N0[channel];
	if (N0[channel] > 0) logprob[channel] += N0[channel] * log(whole_model->getTheta(0));
	logprob[channel] += read_models[channel]->calcLogP();
	for (int i = 0; i < num_threads; ++i) {
		count0[channel] += paramsVecs[channel][i]->count0;
		logprob[channel] += paramsVecs[channel][i]->loglik;
	}
	//  logprob[channel] -= N_eff[channel] * log(whole_model->getProbPass());
	
	if (!keepGoing) whole_model->wrapItUp(count0[channel]);
	else {
		// Run PROBerWholeModel's EM_step procedure
		whole_model->EM_step(count0[channel]);
		
		if (updateReadModel) {
			read_models[channel]->init();
			for (int i = 0; i < num_threads; ++i) read_models[channel]->collect(paramsVecs[channel][i]->estimator);
			read_models[channel]->finish();
		}
	}
}

void EM() {
	int ROUND;
	double prev_logprob, curr_logprob;

	ROUND = 0;
	needCalcConPrb = updateReadModel = true;
	prev_logprob = curr_logprob = -1e300;
	keepGoing = true;

	do {
		++ROUND;

		needCalcConPrb = updateReadModel;
		updateReadModel = needUpdateReadModel(ROUND);

		keepGoing = (ROUND <= MAX_ROUND) && (ROUND <= 2 || (curr_logprob - prev_logprob) / (N_eff[0] + N_eff[1]) > deltaChange);

		// (-) channel
		if (has_control) {
			one_EM_iteration(0, ROUND);
			whole_model->flipState();
		}
		
		// (+) channel
		one_EM_iteration(1, ROUND);
		if (has_control) whole_model->flipState();

		prev_logprob = curr_logprob;
		curr_logprob = logprob[0] + logprob[1];

		if (verbose) printf("Log probability of ROUND %d = %.2f, delta Change = %.10g\n", ROUND - 1, curr_logprob, (curr_logprob - prev_logprob) / (N_eff[0] + N_eff[1]));

	} while (keepGoing);

	if (output_logMAP) {
		char logMAPF[STRLEN];
		sprintf(logMAPF, "%s.logMAP", sampleName);
		FILE *fo = fopen(logMAPF, "w");
		fprintf(fo, "%.0f\n", curr_logprob);
		fclose(fo);
	}
	
	if (verbose) printf("EM is finished!\n");
}

void outputBamFiles(int channel) {
	char inp0F[STRLEN], inpF[STRLEN], inp2F[STRLEN], outF[STRLEN];

	sprintf(outF, "%s_%s.bam", sampleName, channelStr[channel]);
	BamWriter* writer = new BamWriter(outF, hdr, "PROBer");
	AlignmentGroup ag;
	READ_INT_TYPE cnt = 0;
	
	for (int i = 0; i < num_threads; ++i) {
		sprintf(inpF, "%s_%s_%d.bam", imdName, channelStr[channel], i);
		SamParser* parser = new SamParser(inpF, hdr);
		InMemChunk *chunk = paramsVecs[channel][i]->chunk;
		READ_INT_TYPE nreads = chunk->nreads;
		InMemAlignG *a_read = NULL;
		InMemAlign *aligns = NULL;
		
		ag.clear();
		chunk->reset();
		for (READ_INT_TYPE j = 0; j < nreads; ++j) {
			assert(parser->next(ag));
			assert(chunk->next(a_read, aligns));
			
			int size = a_read->size;
			for (int k = 0; k < size; ++k) 
				ag.getAlignment(k)->setFrac(aligns[k].frac);
			writer->write(ag, 2);
			
			++cnt;
			if (verbose && (cnt % 1000000 == 0)) cout<< "Processed "<< cnt<< " reads!"<< endl;
		}
		delete parser;
	}
	
	// write out unalignable reads
	sprintf(inp0F, "%s_%s_N0.bam", imdName, channelStr[channel]);
	SamParser* parser0 = new SamParser(inp0F);
	ag.clear();
	while (parser0->next(ag)) writer->write(ag, 2);
	delete parser0;
	
	// write out filtered reads
	sprintf(inp2F, "%s_%s_N2.bam", imdName, channelStr[channel]);
	SamParser* parser2 = new SamParser(inp2F, hdr);
	ag.clear();
	while (parser2->next(ag)) {
		ag.markAsFiltered(); // Mark each alignment as filtered by append a "ZF:A:!" field
		writer->write(ag, 2);
	}
	delete parser2;
	
	delete writer;
	
	if (verbose) printf("OUTPUT BAM for %s channel is done!\n", channelStr[channel]);
}

void writeResults() {
	// output read model parameters
	char readModelF[STRLEN];

	if (has_control) {
		sprintf(readModelF, "%s_%s.read_model", statName, channelStr[0]);
		read_models[0]->write(readModelF);
	}

	sprintf(readModelF, "%s_%s.read_model", statName, channelStr[1]);
	read_models[1]->write(readModelF);

	
	// output whole model parameters
	whole_model->write(sampleName, statName);

	// output BAM files
	if (output_bam) {
		time_t a = time(NULL);
		// Bam files for (-) channel
		if (has_control) outputBamFiles(0);
		// Bam files for (+) channel
		outputBamFiles(1);

		time_t b = time(NULL);
		char timeF[STRLEN];
		sprintf(timeF, "%s.bam.time", sampleName);
		FILE *fo = fopen(timeF, "w");
		fprintf(fo, "Generating_BAM_files\t%ds\n", int(b - a));
		fclose(fo);
	}

	if (verbose) printf("WriteResults is finished!\n");
}

void release() {
	pthread_attr_destroy(&attr);

	for (int i = 0; i < num_threads; ++i) {
		if (has_control) delete paramsVecs[0][i];
		delete paramsVecs[1][i];
	}

	delete whole_model;
	if (has_control) delete read_models[0];
	delete read_models[1];

	bam_hdr_destroy(hdr);
}

int main(int argc, char* argv[]) {
	if (argc < 7) {
		printf("Usage: PROBer-run-em refName model_type sampleName imdName statName num_of_threads [--read-length read_length] [--maximum-likelihood] [--output-bam] [--output-logMAP] [--no-control] [-q]\n");
		exit(-1);
	}

	strcpy(refName, argv[1]);
	model_type = atoi(argv[2]);
	strcpy(sampleName, argv[3]);
	sprintf(imdName, "%s", argv[4]);
	sprintf(statName, "%s", argv[5]);
	num_threads = atoi(argv[6]);

	output_bam = false;
	output_logMAP = false;
	read_length = -1;
	isMAP = true;
	has_control = true;
	for (int i = 7; i < argc; ++i) {
		if (!strcmp(argv[i], "--read-length")) read_length = atoi(argv[i + 1]);
		if (!strcmp(argv[i], "--maximum-likelihood")) isMAP = false;
		if (!strcmp(argv[i], "--output-bam")) output_bam = true;
		if (!strcmp(argv[i], "--output-logMAP")) output_logMAP = true;
		if (!strcmp(argv[i], "--no-control")) has_control = false;
		if (!strcmp(argv[i], "-q")) verbose = false;
	}

	init();
	EM();
	writeResults();
	release();

	return 0;
}
