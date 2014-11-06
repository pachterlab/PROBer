#include<ctime>
#include<cstdio>
#include<cstring>
#include<cstdlib>
#include<cassert>
#include<fstream>

#include "utils.h"
#include "my_assert.h"
#include "sampling.hpp"

#include "Refs.hpp"

#include "DMSWholeModel.hpp"
#include "DMSReadModel.hpp"

using namespace std;

Refs refs;

DMSWholeModel *whole_model;
DMSReadModel *read_model;

seedType seed;
Sampler *sampler;

int N;
int model_type;

char refF[STRLEN], outF1[STRLEN], outF2[STRLEN];
ofstream out1, out2;

int main(int argc, char* argv[]) {
  if (argc < 8 || argc > 9) {
    printf("Usage: dms-seq-simulate-reads reference_name config_file whole_model_input_name read_model_file <'minus' or 'plus'> number_of_reads output_name [seed]\n");
    return 0;
  }

  N = atoi(argv[6]);
  if (argc == 9) { 
    int len = strlen(argv[8]);
    seed = 0;
    for (int i = 0; i < len; ++i) seed = seed * 10 + (argv[8][i] - '0');
  }
  else seed = time(NULL);

  sampler = new Sampler(seed);
  
  whole_model = new DMSWholeModel(argv[2]);
  whole_model->read(argv[3]);
  if (!strcmp(argv[5], "plus")) whole_model->read(argv[3]);

  sprintf(refF, "%s.seq", argv[1]);
  refs.loadRefs(refF);

  read_model = new DMSReadModel(&refs, sampler);
  read_model->read(argv[4]);

  model_type = read_model->getModelType();

  if (model_type < 2) {
    sprintf(outF1, "%s_%s.fq", argv[7], argv[5]);
    out1.open(outF1);
  }
  else {
    sprintf(outF1, "%s_%s_1.fq", argv[7], argv[5]);
    sprintf(outF2, "%s_%s_2.fq", argv[7], argv[5]);
    out1.open(outF1);
    out2.open(outF2);
  }

  int tid, pos, frag_len;
  whole_model->startSimulation();
  read_model->startSimulation();
  for (int i = 0; i < N; ++i) {
    whole_model->simulate(sampler, tid, pos, frag_len);
    if (model_type < 2) read_model->simulate(i, tid, pos, frag_len, &out1);
    else read_model->simulate(i, tid, pos, frag_len, &out1, &out2);
    if ((i + 1) % 1000000 == 0) printf("GEN %d!\n", i + 1);
  }
  whole_model->finishSimulation();
  read_model->finishSimulation();

  out1.close();
  if (model_type >= 2) out2.close();

  delete whole_model;
  delete read_model;
  delete sampler;
  
  return 0;
}
