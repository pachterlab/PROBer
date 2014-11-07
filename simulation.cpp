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
    printf("Usage: dms-seq-simulate-reads reference_name config_file whole_model_input_name read_model_file channel('minus' or 'plus') number_of_reads output_name [seed]\n\n");
    printf("Description:\n");
    printf("  This program simulate reads using parameters learned from data by 'dms-seq-estimate-parameters'.\n\n");
    printf("Arguments:\n");
    printf("  reference_name: The reference's name, should be same as the ones used in 'dms-seq-prepare-reference' and 'dms-seq-estimate-parameters'.\n");
    printf("  config_file: Contains primer length, size selection min and max fragment size etc. 'sample_name.temp/sample_name_minus.config' and 'sample_name.temp/sample_name_plus.config' can be used here.\n");
    printf("  whole_model_input_name: This should be the 'sample_name' used in 'dms-seq-estimate-parameters'.\n");
    printf("  read_mode_file: 'sample_name.stat/sample_name_minus.read_model' if chanel is 'minus', 'sample_name.stat/sample_name_plus.read_model' if channel is 'plus'.\n");
    printf("  channel: 'minus' for minus channel and 'plus' for plus channel, no quotation.\n");
    printf("  number_of_reads: Number of reads to simulate.\n");
    printf("  output_name: Output file prefix. Simulated single-end reads will be written into 'output_name_[minus/plus].[fa/fq]'. Simulated paired-end reads will be written into 'output_name_[minus/plus]_1.[fa/fq]' and 'output_name_[minus/plus]_2.[fa/fq].\n");
    printf("  [seed]: Optional argument, the seed used for simulation.\n");
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
    sprintf(outF1, "%s_%s.%s", argv[7], argv[5], (model_type == 0 ? "fa" : "fq"));
    out1.open(outF1);
  }
  else {
    sprintf(outF1, "%s_%s_1.%s", argv[7], argv[5], (model_type == 2 ? "fa" : "fq"));
    sprintf(outF2, "%s_%s_2.%s", argv[7], argv[5], (model_type == 2 ? "fa" : "fq"));
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
