#include<cmath>
#include<cstdio>
#include<cstring>
#include<cstdlib>
#include<cassert>
#include<vector>
#include<string>
#include<iostream>
#include<pthread.h>

#include "utils.h"
#include "my_assert.h"

#include "Refs.hpp"
#include "Transcripts.hpp"

#include "BamAlignment.hpp"
#include "AlignmentGroup.hpp"
#include "SamParser.hpp"
#include "BamWriter.hpp"

#include "DMSWholeModel.hpp"
#include "DMSReadModel.hpp"
#include "InMemoryStructs.hpp"

using namespace std;

const int MAX_ROUND = 200;

// Parameter struct to pass parameters to each subprocess
struct InMemParams {
  int no; // thread number

  DMSWholeModel *whole_model;
  DMSReadModel *read_model;

  DMSReadModel *estimator; // slave model that is used to estimate model parameters
  InMemChunk *chunk; // A chunk of memory to record all in-memory information for reads and alignments associated with this thread

  double count0; // sum of noise read fractions
  double loglik; // log likelihood

  InMemParams(int no, DMSWholeModel* whole_model, DMSReadModel* read_model, READ_INT_TYPE nreads, HIT_INT_TYPE nlines) {
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

int M; // Number of transcripts
int N0, N_eff; // Number of unalignable reads, number of effective reads (unaligned + aligned)
int model_type; 
int num_threads;
int read_length;
bool isMAP;

char refName[STRLEN], sampleName[STRLEN], imdName[STRLEN], statName[STRLEN], channel[STRLEN];

Refs refs;
Transcripts transcripts;

DMSWholeModel *whole_model;
DMSReadModel *read_model;

bool needCalcConPrb, updateReadModel;
vector<InMemParams*> paramsVec;
vector<pthread_t> threads;
pthread_attr_t attr;
int rc;

bool output_bam;

void init() {
  char refF[STRLEN], tiF[STRLEN];
  char configF[STRLEN];
  char bamF[STRLEN], partitionF[STRLEN];

  // Load references
  sprintf(refF, "%s.seq", refName);
  refs.loadRefs(refF);
  M = refs.getM();
  
  sprintf(tiF, "%s.ti", refName);
  transcripts.readFrom(tiF);
  transcripts.buildMappings(imdName);

  // Create DMSWholeModel
  sprintf(configF, "%s.config", imdName);
  whole_model = new DMSWholeModel(configF, &transcripts, num_threads, read_length, isMAP);
  if (!strcmp(channel, "plus")) whole_model->read(sampleName); // Read gamma because this run is used to estimate betas

  // Create DMSReadModel
  read_model = new DMSReadModel(model_type, &refs, read_length);

  // Preprocess reads and alignments
  SamParser *parser = NULL;
  AlignmentGroup ag;

  // Preprocess unalignable reads
  sprintf(bamF, "%s_N0.bam", imdName);
  parser = new SamParser('b', bamF, NULL);

  N0 = 0;
  while (parser->next(ag)) {
    read_model->update_preprocess(ag, false);
    ++N0;
  }
  delete parser;
  if (verbose) { printf("Unalignable reads are preprocessed!\n"); }

  // Preprocess alignable reads
  sprintf(partitionF, "%s.partition", imdName);
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

  N_eff = N0;
  paramsVec.assign(num_threads, NULL);
  for (int i = 0; i < num_threads; ++i) {
    fin>> id>> nreads>> nlines;
    assert(id == i);
    N_eff += nreads;
    paramsVec[i] = new InMemParams(i, whole_model, read_model, nreads, nlines);

    sprintf(bamF, "%s_%d.bam", imdName, i);
    parser = new SamParser('b', bamF, NULL);
    rid = 0;
    ag.clear();
    while (parser->next(ag)) {
      is_paired = ag.isPaired();
      seqlen = !is_paired ? ag.getSeqLength() : 0; 

      assert(paramsVec[i]->chunk->next(a_read, aligns));
      a_read->size = ag.size();
      
      for (int j = 0; j < ag.size(); ++j) {
	BamAlignment *ba = ag.getAlignment(j);
	aligns[j].tid = transcripts.getInternalSid(ba->getTid());
	aligns[j].pos = ba->getLeftMostPos();

	if (is_paired) aligns[j].fragment_length = ba->getInsertSize();
	else if (seqlen < read_length) aligns[j].fragment_length = seqlen;
	else aligns[j].fragment_length = 0;

	aligns[j].frac = 1.0 / ag.size();
      }
      whole_model->addAlignments(a_read, aligns);
      read_model->update_preprocess(ag, true);
      ++rid;

      if (verbose && (rid % 1000000 == 0)) cout<< "Loaded "<< rid<< " reads!"<< endl;
    }
    assert(rid == nreads);

    delete parser;
    if (verbose) { printf("Thread %d's data is preprocessed!\n", i); }

    chunk = paramsVec[i]->chunk;
    for (HIT_INT_TYPE j = 0; j < nlines; ++j) 
      if (chunk->aligns[j].conprb == -1.0) ++cnt;
  }

  if (verbose) { printf("There are %d alignments filtered!\n", cnt); }

  for (int i = 0; i < num_threads; ++i) 
    paramsVec[i]->estimator = new DMSReadModel(read_model);
  whole_model->init_for_EM();
  read_model->finish_preprocess();

  threads.assign(num_threads, pthread_t());
  /* set thread attribute to be joinable */
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

  if (verbose) { printf("Preprocess data is finished!\n"); }
}

void* E_STEP(void* arg) {
  InMemParams *params = (InMemParams*)arg;
  DMSWholeModel *whole_model = params->whole_model;
  DMSReadModel *read_model = params->read_model;
  DMSReadModel *estimator = params->estimator;
  InMemChunk *chunk = params->chunk;

  SamParser *parser = NULL;
  AlignmentGroup ag;

  params->count0 = 0.0;
  params->loglik = 0.0;

  chunk->reset();

  if (needCalcConPrb || updateReadModel) {
    char bamF[STRLEN];
    sprintf(bamF, "%s_%d.bam", imdName, params->no);
    parser = new SamParser('b', bamF, NULL); 
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

void EM() {
  int ROUND;
  double count0;
  double loglik;

  ROUND = 0;
  needCalcConPrb = updateReadModel = true;
  loglik = 0.0;

  do {
    ++ROUND;

    needCalcConPrb = updateReadModel;
    updateReadModel = needUpdateReadModel(ROUND);
    
    // E step
    for (int i = 0; i < num_threads; ++i) {
      rc = pthread_create(&threads[i], &attr, E_STEP, (void*)paramsVec[i]);
      pthread_assert(rc, "pthread_create", "Cannot create thread " + itos(i) + " (numbered from 0) at EM ROUND " + itos(ROUND) + "!");
    }
    
    for (int i = 0; i < num_threads; ++i) {
      rc = pthread_join(threads[i], NULL);
      pthread_assert(rc, "pthread_join", "Cannot join thread " + itos(i) + " (numbered from 0) at EM ROUND " + itos(ROUND) + "!");
    }

    count0 = N0;
    loglik = N0 * log(whole_model->getTheta(0)) + read_model->calcLogP();
    for (int i = 0; i < num_threads; ++i) {
      count0 += paramsVec[i]->count0;
      loglik += paramsVec[i]->loglik;
    }
    loglik -= N_eff * log(whole_model->getProbPass());

    if (verbose) printf("E step is done. Loglik of ROUND %d is: %.2f\n", ROUND - 1, loglik);

    if (ROUND > MAX_ROUND) {
      whole_model->update(count0);
      continue;
    }

    // Run DMSWholeModel's runEM procedure
    whole_model->runEM(count0);

    if (updateReadModel) {
      read_model->init();
      for (int i = 0; i < num_threads; ++i) read_model->collect(paramsVec[i]->estimator);
      read_model->finish();
    }

    if (verbose) printf("ROUND %d finished!\n", ROUND);

  } while (ROUND <= MAX_ROUND);
  
  if (verbose) printf("EM is finished!\n");
}

void writeResults() {
  // output read model parameters
  char readModelF[STRLEN];
  sprintf(readModelF, "%s.read_model", statName);
  read_model->write(readModelF);
  
  // output whole model parameters
  whole_model->setDefault();
  whole_model->write(sampleName);

  // output BAM files
  if (output_bam) {
    char inp0F[STRLEN], inpF[STRLEN], inp2F[STRLEN], outF[STRLEN];
    sprintf(inp0F, "%s_N0.bam", imdName);
    sprintf(outF, "%s_%s.transcripts.bam", sampleName, channel);
    
    SamParser* parser0 = new SamParser('b', inp0F, NULL);
    BamWriter* writer = new BamWriter(outF, parser0->getHeader(), "DMS-Seq");
    AlignmentGroup ag;
    READ_INT_TYPE cnt = 0;

    for (int i = 0; i < num_threads; ++i) {
      sprintf(inpF, "%s_%d.bam", imdName, i);
      SamParser* parser = new SamParser('b', inpF, NULL);
      InMemChunk *chunk = paramsVec[i]->chunk;
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
    ag.clear();
    while (parser0->next(ag)) writer->write(ag, 2);
    delete parser0;

    // write out filtered reads
    sprintf(inp2F, "%s_N2.bam", imdName);
    SamParser* parser2 = new SamParser('b', inp2F, NULL);
    ag.clear();
    while (parser2->next(ag)) {
      ag.markAsFiltered(); // Mark each alignment as filtered by append a "ZF:A:!" field
      writer->write(ag, 2);
    }
    delete parser2;

    delete writer;
    
    if (verbose) printf("OUTPUT BAM is written!\n");
  }

  if (verbose) printf("WriteResults is finished!\n");
}

void release() {
  pthread_attr_destroy(&attr);

  for (int i = 0; i < num_threads; ++i) delete paramsVec[i];
  delete whole_model;
  delete read_model;
}

int main(int argc, char* argv[]) {
  if (argc < 8) {
    printf("Usage: dms-seq-run-em refName model_type sampleName imdName statName channel<'minus' or 'plus'> num_of_threads [--read-length read_length] [--MAP] [--output-bam] [-q]\n");
    exit(-1);
  }

  strcpy(refName, argv[1]);
  model_type = atoi(argv[2]);
  strcpy(sampleName, argv[3]);
  sprintf(imdName, "%s_%s", argv[4], argv[6]);
  sprintf(statName, "%s_%s", argv[5], argv[6]);
  strcpy(channel, argv[6]);
  num_threads = atoi(argv[7]);

  verbose = true;
  output_bam = false;
  read_length = -1;
  isMAP = false;
  for (int i = 8; i < argc; ++i) {
    if (!strcmp(argv[i], "--read-length")) read_length = atoi(argv[i + 1]);
    if (!strcmp(argv[i], "--MAP")) isMAP = true;
    if (!strcmp(argv[i], "--output-bam")) output_bam = true;
    if (!strcmp(argv[i], "-q")) verbose = false;
  }

  init();
  EM();
  writeResults();
  release();

  return 0;
}
