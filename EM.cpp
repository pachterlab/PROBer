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
};

int M; // Number of transcripts
char refF[STRLEN], tiF[STRLEN];

DMSWholeModel *whole_model;
DMSReadModel *read_model;

vector<InMemParams> paramsVec;
vector<pthread_t> threads;
pthread_attr_t attr;
int rc;

int main(int argc, char* argv[]) {
  if (argc < 6) {
    printf("Usage: dms-seq-run-em refName sampleName imdName statName num_of_threads [-q]\n");
  }

  verbose = true;
  for (int i = 6; i < argc; ++i) {
    if (!strcmp(argv[i], "-q")) verbose = false;
  }

  sprintf(refF, "%s.seq", argv[1]);
  refs.loadRefs(refF);
  M = refs.getM();

  sprintf(tiF, "%s.ti", argv[1]);
  transcripts.readFrom(tiF);


  return 0;
}
