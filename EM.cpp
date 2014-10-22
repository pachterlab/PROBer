#include<cmath>
#include<cstdio>
#include<cstring>
#include<cstdlib>
#include<cassert>
#include<vector>
#include<string>
#include<pthread.h>

#include "utils.h"
#include "my_assert.h"

#include "Transcripts.hpp"
#include "Refs.hpp"

#include "AlignmentGroup.hpp"
#include "SamParser.hpp"

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
  vector<InMemAlignG> ags; // alignment groups

  vector<double> fracs; // fractions for each alignment in a read
  double count0; // sum of noise read fractions
  double loglik; // log likelihood

  InMemParams(int no, DMSWholeModel* whole_model, DMSReadModel* read_model, READ_INT_TYPE nreads) {
    this->no = no;
    this->whole_model = whole_model;
    this->read_model = read_model;
    estimator = NULL;
    ags.assign(nreads, InMemAlignG());
    fracs.clear();
    count0 = loglik = 0.0;
  }

  ~InMemParams() {
    ags.clear();
    if (estimator != NULL) delete estimator;
  }
};

int M; // Number of transcripts
int N0, N_eff; // Number of unalignable reads, number of effective reads (unaligned + aligned)
int model_type; 
int num_threads;

char refName[STRLEN], sampleName[STRLEN], imdName[STRLEN], statName[STRLEN];

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
  char cntF[STRLEN];

  // Load references
  sprintf(refF, "%s.seq", refName);
  refs.loadRefs(refF);
  M = refs.getM();
  
  sprintf(tiF, "%s.ti", refName);
  transcripts.readFrom(tiF);
  transcripts.buildMappings(imdName);

  // Create DMSWholeModel
  sprintf(configF, "%s.config", imdName);
  whole_model = new DMSWholeModel(configF, &transcripts, num_threads);

  // Create DMSReadModel
  read_model = new DMSReadModel(model_type, refs);

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
  int max_size; // maximum alignment sizes in a read for a thread

  N_eff = N0;
  paramsVec.assign(num_threads, NULL);
  for (int i = 0; i < num_threads; ++i) {
    fin>> id>> nreads>> nlines;
    N_eff += nreads;
    paramsVec[i] = new InMemParams(i, whole_model, read_model, nreads);

    sprintf(bamF, "%s_%d.bam", imdName, i);
    parser = new SamParser('b', bamF, NULL);
    rid = 0;
    max_size = 0;
    while (parser->next(ag)) {
      paramsVec[i]->ags[rid].size = ag.size();
      if (max_size < ag.size()) max_size = ag.size();
      paramsVec[i]->ags[rid].aligns = new InMemAlign*[ag.size()];
      for (int j = 0; j < ag.size(); ++j) {
	BamAlignment *ba = ag.getAlignment(j);
	paramsVec[i]->ags[rid].aligns[j] = new InMemAlign(transcripts.getInternalSid(ba->getTid()), ba->getLeftMostPos(), (ba->isPaired() ? ba->getInsertSize() : 0), 1.0 / ag.size());
      }
      whole_model->addAlignments(paramsVec[i]->ags[rid]);
      read_model->update_preprocess(ag, true);
      ++rid;
    }
    delete parser;
    paramsVec[i]->fracs.assign(max_size + 1, 0.0);
    if (verbose) { printf("Thread %d's data is preprocessed!\n", i); }
  }

  for (int i = 0; i < num_threads; ++i) 
    paramsVec[i]->estimator = new DMSReadModel(read_model);

  whole_model->init_for_EM();

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
  vector<InMemAlignG>& ags = params->ags;
  vector<double>& fracs = params->fracs;

  SamParser *parser = NULL;
  AlignmentGroup ag;

  params->count0 = 0.0;
  params->loglik = 0.0;

  if (needCalcConPrb || updateReadModel) {
    char bamF[STRLEN];
    sprintf(bamF, "%s_%d.bam", imdName, params->no);
    parser = new SamParser('b', bamF, NULL); 
  }
  if (updateReadModel) estimator->init();

  READ_INT_TYPE nreads = ags.size();
  int size;
  double sum;

  for (READ_INT_TYPE i = 0; i < nreads; ++i) {
    if (needCalcConPrb || updateReadModel) {
      assert(parser->next(ag));
    }
    
    if (needCalcConPrb) read_model->setProbs(ags[i], ag);

    size = ags[i].size;
    sum = 0.0;
    for (int j = 0; j < size; ++j) {
      fracs[j] = whole_model->getProb(ags[i].aligns[j]->tid, ags[i].aligns[j]->pos, ags[i].aligns[j]->fragment_length) * ags[i].aligns[j]->frac;
      sum += fracs[j];
    }
    fracs[size] = whole_model->getProb(0) * ags[i].noise_prob;
    sum += fracs[size];
    assert(sum > 0.0);
    params->loglik += log(sum);
    for (int j = 0; j <= size; ++j) fracs[j] /= sum;
    params->count0 += fracs[size];

    if (updateReadModel) estimator->update(ags[i], ag, fracs);    
  }

  if (parser != NULL) delete parser;
}

inline bool needUpdateReadModel(int ROUND) {
  return ROUND <= 10;
}

void EM() {
  int ROUND;
  double count0;
  double loglik, loglik_old;

  ROUND = 0;
  needCalcConPrb = updateReadModel = true;

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
    loglik_old = loglik;
    loglik = N0 * log(whole_model->getTheta(0)) + read_model->calcLogP();
    for (int i = 0; i < num_threads; ++i) {
      count0 += paramsVec[i]->count0;
      loglik += paramsVec[i]->loglik;
    }
    loglik -= N_eff * log(whole_model->getProbPass());

    printf("Main EM finished. Loglik of current parameters is: %.2f\n", loglik);

    // Run DMSWholeModel's runEM procedure
    whole_model->runEM(count0, 1);

    if (updateReadModel) {
      read_model->init();
      for (int i = 0; i < num_threads; ++i) read_model->collect(paramsVec[i]->estimator);
      read_model->finish();
    }

  } while (ROUND < MAX_ROUND);
}

void* calc_expected_alignments(void* arg) {
  InMemParams *params = (InMemParams*)arg;
  DMSWholeModel *whole_model = params->whole_model;
  DMSReadModel *read_model = params->read_model;
  vector<InMemAlignG>& ags = params->ags;
  vector<double>& fracs = params->fracs;

  READ_INT_TYPE nreads = ags.size();
  int size;
  double sum;

  for (READ_INT_TYPE i = 0; i < nreads; ++i) {
    size = ags[i].size;
    sum = 0.0;
    for (int j = 0; j < size; ++j) {
      fracs[j] = whole_model->getProb(ags[i].aligns[j]->tid, ags[i].aligns[j]->pos, ags[i].aligns[j]->fragment_length) * ags[i].aligns[j]->frac;
      sum += fracs[j];
    }
    fracs[size] = whole_model->getProb(0) * ags[i].noise_prob;
    sum += fracs[size];
    assert(sum > 0.0);
    for (int j = 0; j < size; ++j) ags[i].aligns[j]->frac = fracs[j] / sum;
    ags[i].noise_prob = fracs[size] / sum;
  }
}

void writeResults() {
  if (output_bam) {
    // Calculate expected alignment weights
    for (int i = 0; i < num_threads; ++i) {
      rc = pthread_create(&threads[i], &attr, calc_expected_alignments, (void*)paramsVec[i]);
      pthread_assert(rc, "pthread_create", "Cannot create thread " + itos(i) + " (numbered from 0) at EM ROUND " + itos(ROUND) + "!");
    }
    
    for (int i = 0; i < num_threads; ++i) {
      rc = pthread_join(threads[i], NULL);
      pthread_assert(rc, "pthread_join", "Cannot join thread " + itos(i) + " (numbered from 0) at EM ROUND " + itos(ROUND) + "!");
    }

    char inp0F[STRLEN], inpF[STRLEN], outF[STRLEN];
    sprintf(inp0F, "%s_N0.bam", imdName);
    sprintf(outF, "%s.transcripts.bam", sampleName);
    
    SamParser* parser0 = new SamParser('b', inp0F, NULL);
    BamWriter* writer = new BamWriter(outF, parser0->getHeader(), "DMS-Seq");
    AlignmentGroup ag;
    for (int i = 0; i < num_threads; ++i) {
      sprintf(inpF, "%s_%d.bam", imdName, i);
      SamParser* parser = new SamParser('b', inpF, NULL);
      vector<InMemAlignG>& ags = paramsVec[i]->ags;
      READ_INT_TYPE nreads = ags.size();
      for (READ_INT_TYPE j = 0; j < nreads; ++j) {
	assert(parser->next(ag));
	int size = ags[j].size;
	for (int k = 0; k < size; ++k) 
	  ag.getAlignment(k)->setFrac(ags[j].aligns[k]->frac);
	writer->write(ag, 2);
      }
      delete parser;
    }
    // write out unalignable reads
    while (parser0->next(ag)) writer->write(ag, 2);

    delete parser0;
    delete writer;
  }

  // output read model parameters
  char readModelF[STRLEN];
  sprintf(readModelF, "%s.read_model", statName);
  read_model->write(readModelF);
  
  // output whole model parameters
  whole_model->write(sampleName);
}

void release() {
  pthread_attr_destroy(&attr);

  for (int i = 0; i < num_threads; ++i) delete paramsVec[i];
  delete whole_model;
  delete read_model;
}

int main(int argc, char* argv[]) {
  if (argc < 7) {
    printf("Usage: dms-seq-run-em refName model_type sampleName imdName statName num_of_threads [--output-bam] [-q]\n");
  }

  strcpy(refName, argv[1]);
  model_type = atoi(argv[2]);
  strcpy(sampleName, argv[3]);
  strcpy(imdName, argv[4]);
  strcpy(statName, argv[5]);
  num_threads = atoi(argv[6]);

  verbose = true;
  output_bam = false;
  for (int i = 7; i < argc; ++i) {
    if (!strcmp(argv[i], "--output-bam")) output_bam = true;
    if (!strcmp(argv[i], "-q")) verbose = false;
  }

  init();
  EM();
  writeResults();
  release();

  return 0;
}
