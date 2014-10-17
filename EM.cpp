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

// Parameter struct to pass parameters to each subprocess
struct InMemParams {
  DMSWholeModel *whole_model;
  DMSReadModel *read_model;

  DMSReadModel *estimator; // slave model that is used to estimate model parameters
  std::vector<InMemAlignG> ags; // alignment groups

  InMemParams(DMSWholeModel* whole_model, DMSReadModel* read_model, READ_INT_TYPE nreads) {
    this->whole_model = whole_model;
    this->read_model = read_model;
    estimator = NULL;
    ags.assign(nreads, InMemAlignG());
  }

  ~InMemParams() {
    if (estimator != NULL) delete estimator;
  }
};

int M; // Number of transcripts
int model_type; 
int num_threads;

char refName[STRLEN], sampleName[STRLEN], imdName[STRLEN], statName[STRLEN];

DMSWholeModel *whole_model;
DMSReadModel *read_model;

vector<InMemParams*> paramsVec;
vector<pthread_t> threads;
pthread_attr_t attr;
int rc;

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
  whole_model = new DMSWholeModel(configF, &transcripts, num_threads);

  // Create DMSReadModel
  read_model = new DMSReadModel(model_type, refs);

  // Preprocess reads and alignments
  SamParser *parser = NULL;
  AlignmentGroup ag;

  // Preprocess unalignable reads
  sprintf(bamF, "%s_N0.bam", imdName);
  parser = new SamParser('b', bamF, NULL);

  while (parser->next(ag)) {
    read_model->update_preprocess(ag, false);
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

  paramsVec.assign(num_threads, NULL);
  for (int i = 0; i < num_threads; ++i) {
    fin>> id>> nreads>> nlines;
    paramsVec[i] = new InMemParams(whole_model, read_model, nreads);

    sprintf(bamF, "%s_%d.bam", imdName, i);
    parser = new SamParser('b', bamF, NULL);
    rid = 0;
    while (parser->next(ag)) {
      paramsVec[i]->ags[rid].size = ag.size();
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
    if (verbose) { printf("Thread %d's data is preprocessed!\n", i); }
  }

  for (int i = 0; i < num_threads; ++i) 
    paramsVec[i]->estimator = new DMSReadModel(read_model);

  whole_model->init_for_EM();

  if (verbose) { printf("Preprocess data is finished!\n"); }
}

int main(int argc, char* argv[]) {
  if (argc < 7) {
    printf("Usage: dms-seq-run-em refName model_type sampleName imdName statName num_of_threads [-q]\n");
  }

  strcpy(refName, argv[1]);
  model_type = atoi(argv[2]);
  strcpy(sampleName, argv[3]);
  strcpy(imdName, argv[4]);
  strcpy(statName, argv[5]);
  num_threads = atoi(argv[6]);

  verbose = true;
  for (int i = 7; i < argc; ++i) {
    if (!strcmp(argv[i], "-q")) verbose = false;
  }

  init();

  return 0;
}
