#include<cmath>
#include<ctime>
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

#include "PROBerWholeModel.hpp"
#include "PROBerReadModel.hpp"
#include "InMemoryStructs.hpp"

using namespace std;

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

int MAX_ROUND;

int M; // Number of transcripts
int N0[2], N_eff[2]; // Number of unalignable reads, number of effective reads (unaligned + aligned)
double count0[2], loglik[2]; // Used in EM algorithm, number of unalignable reads and log likelihood for each channel

int model_type; 
int num_threads;
int read_length;
bool isMAP;

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

bool output_bam;

// Preprocess reads and alignments
void preprocessAlignments(int channel) {
  char bamF[STRLEN], partitionF[STRLEN];
  SamParser *parser = NULL;
  AlignmentGroup ag;

  if (verbose) { printf("Begin to preprocess BAM files for channel %s!\n", channelStr[channel]); }

  // Preprocess unalignable reads
  sprintf(bamF, "%s_%s_N0.bam", imdName, channelStr[channel]);
  parser = new SamParser('b', bamF, NULL);

  N0[channel] = 0;
  while (parser->next(ag)) {
    read_models[channel]->update_preprocess(ag, false);
    ++N0[channel];
  }
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
    parser = new SamParser('b', bamF, NULL);
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
	else if (seqlen < read_length) aligns[j].fragment_length = seqlen;
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
  whole_model->flipState();
}

void init() {
  char refF[STRLEN], tiF[STRLEN];
  char configF[STRLEN];

  // Load references
  sprintf(refF, "%s.seq", refName);
  refs.loadRefs(refF);
  M = refs.getM();
  
  sprintf(tiF, "%s.ti", refName);
  transcripts.readFrom(tiF);
  char imd_name[STRLEN];
  sprintf(imd_name, "%s_minus", imdName);
  transcripts.buildMappings(imd_name);

  // Create PROBerWholeModel
  sprintf(configF, "%s.config", imdName);
  whole_model = new PROBerWholeModel(configF, 2, &transcripts, num_threads, read_length, isMAP);

  // Create PROBerReadModels
  read_models[0] = new PROBerReadModel(model_type, &refs, read_length);
  read_models[1] = new PROBerReadModel(model_type, &refs, read_length);

  // Preprocess data for (-)
  preprocessAlignments(0);
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
    sprintf(bamF, "%s_%s_%d.bam", imdName, channelStr[whole_model->getChannel()], params->no);
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

void one_EM_iteration(int channel, int ROUND) {
  assert(whole_model->getChannel() == channel);

  // init
  if (ROUND == 1) whole_model->init();

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
  loglik[channel] = N0[channel] * log(whole_model->getTheta(0)) + read_models[channel]->calcLogP();
  for (int i = 0; i < num_threads; ++i) {
    count0[channel] += paramsVecs[channel][i]->count0;
    loglik[channel] += paramsVecs[channel][i]->loglik;
  }
  loglik[channel] -= N_eff[channel] * log(whole_model->getProbPass());
  
  if (ROUND > MAX_ROUND) whole_model->wrapItUp(count0[channel]);
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
  double prev_loglik, curr_loglik;

  ROUND = 0;
  needCalcConPrb = updateReadModel = true;
  prev_loglik = curr_loglik = -1e300;

  do {
    ++ROUND;

    needCalcConPrb = updateReadModel;
    updateReadModel = needUpdateReadModel(ROUND);

    // (-) channel
    one_EM_iteration(0, ROUND);
    whole_model->flipState();

    // (+) channel
    one_EM_iteration(1, ROUND);
    whole_model->flipState();

    prev_loglik = curr_loglik;
    curr_loglik = loglik[0] + loglik[1];

    //if (verbose) printf("Loglik of ROUND %d is: %.2f\n", ROUND - 1, loglik[0] + loglik[1]);
    if (verbose) printf("Loglik of ROUND %d is: %.2f\n", ROUND - 1, curr_loglik);
    if (ROUND > 1) printf("delta_loglik = %.10g, avg_delta_loglik = %.10g\n", (curr_loglik - prev_loglik), (curr_loglik - prev_loglik) / (N_eff[0] + N_eff[1]));

  } while (ROUND <= MAX_ROUND);
  
  if (verbose) printf("EM is finished!\n");
}

void outputBamFiles(int channel) {
  char inp0F[STRLEN], inpF[STRLEN], inp2F[STRLEN], outF[STRLEN];

  sprintf(inp0F, "%s_%s_N0.bam", imdName, channelStr[channel]);
  sprintf(outF, "%s_%s.transcripts.bam", sampleName, channelStr[channel]);
  
  SamParser* parser0 = new SamParser('b', inp0F, NULL);
  BamWriter* writer = new BamWriter(outF, parser0->getHeader(), "PROBer");
  AlignmentGroup ag;
  READ_INT_TYPE cnt = 0;
  
  for (int i = 0; i < num_threads; ++i) {
    sprintf(inpF, "%s_%s_%d.bam", imdName, channelStr[channel], i);
    SamParser* parser = new SamParser('b', inpF, NULL);
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
  ag.clear();
  while (parser0->next(ag)) writer->write(ag, 2);
  delete parser0;
  
  // write out filtered reads
  sprintf(inp2F, "%s_%s_N2.bam", imdName, channelStr[channel]);
  SamParser* parser2 = new SamParser('b', inp2F, NULL);
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
  for (int i = 0; i < 2; ++i) {
    sprintf(readModelF, "%s_%s.read_model", statName, channelStr[i]);
    read_models[i]->write(readModelF);
  }
  
  // output whole model parameters
  whole_model->write(sampleName, statName);

  // output BAM files
  if (output_bam) {
    time_t a = time(NULL);
    // Bam files for (-) channel
    outputBamFiles(0);
    // Bam files for (+) channel
    outputBamFiles(1);
    time_t b = time(NULL);

    char timeF[STRLEN];
    sprintf(timeF, "%s.bam.time", sampleName);
    FILE *fo = fopen(timeF, "w");
    fprintf(fo, "%ds or %.2fm or %.2fh.\n", int(b - a), (b - a) / 60.0, (b - a) / 3600.0);
    fclose(fo);
  }

  if (verbose) printf("WriteResults is finished!\n");
}

void release() {
  pthread_attr_destroy(&attr);

  for (int i = 0; i < num_threads; ++i) {
    delete paramsVecs[0][i];
    delete paramsVecs[1][i];
  }

  delete whole_model;
  delete read_models[0];
  delete read_models[1];
}

int main(int argc, char* argv[]) {
  if (argc < 7) {
    printf("Usage: PROBer-run-em refName model_type sampleName imdName statName num_of_threads [--read-length read_length] [--maximum-likelihood] [--output-bam] [--rounds rounds] [-q]\n");
    exit(-1);
  }

  strcpy(refName, argv[1]);
  model_type = atoi(argv[2]);
  strcpy(sampleName, argv[3]);
  sprintf(imdName, "%s", argv[4]);
  sprintf(statName, "%s", argv[5]);
  num_threads = atoi(argv[6]);

  verbose = true;
  output_bam = false;
  read_length = -1;
  isMAP = true;
  MAX_ROUND = 400;
  for (int i = 7; i < argc; ++i) {
    if (!strcmp(argv[i], "--read-length")) read_length = atoi(argv[i + 1]);
    if (!strcmp(argv[i], "--maximum-likelihood")) isMAP = false;
    if (!strcmp(argv[i], "--output-bam")) output_bam = true;
    if (!strcmp(argv[i], "--rounds")) MAX_ROUND = atoi(argv[i + 1]);
    if (!strcmp(argv[i], "-q")) verbose = false;
  }

  init();
  EM();
  writeResults();
  release();

  return 0;
}
