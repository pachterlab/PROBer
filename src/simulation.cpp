#include<ctime>
#include<cstdio>
#include<cstring>
#include<cstdlib>
#include<cassert>
#include<fstream>
#include<string>

#include "utils.h"
#include "my_assert.h"
#include "sampling.hpp"

#include "Refs.hpp"

#include "PROBerWholeModel.hpp"
#include "PROBerReadModel.hpp"

using namespace std;

bool verbose = true; // define verbose


Refs refs;

PROBerWholeModel *whole_model;
PROBerReadModel *read_model;

seedType seed;
Sampler *sampler;

int M, N;
int model_type;

int sim_tid; // tid used for simulation

char statName[STRLEN];
char refF[STRLEN], readModelF[STRLEN], outF1[STRLEN], outF2[STRLEN];
ofstream out1, out2;

bool has_control;

int main(int argc, char* argv[]) {
  if (argc < 7) {
    printf("Usage: PROBer-simulate-reads reference_name config_file sample_name channel('minus' or 'plus') number_of_reads output_name [--seed seed] [--transcript name] [--read-model-file read_model_file] [--no-control]\n\n");
    printf("Description:\n");
    printf("  This program simulate reads using parameters learned from data by 'PROBer-estimate-parameters'.\n\n");
    printf("Arguments:\n");
    printf("  reference_name: The reference's name, should be same as the ones used in 'PROBer-prepare-reference' and 'PROBer-estimate-parameters'.\n");
    printf("  config_file: Contains primer length, size selection min and max fragment size etc. 'sample_name.temp/sample_name_minus.config' and 'sample_name.temp/sample_name_plus.config' can be used here.\n");
    printf("  sample_name: This should be the 'sample_name' used in 'PROBer-estimate-parameters'. No slash should be in the end of this string.\n");
    printf("  channel: 'minus' for minus channel and 'plus' for plus channel, no quotation.\n");
    printf("  number_of_reads: Number of reads to simulate.\n");
    printf("  output_name: Output file prefix. Simulated single-end reads will be written into 'output_name_[minus/plus].[fa/fq]'. Simulated paired-end reads will be written into 'output_name_[minus/plus]_1.[fa/fq]' and 'output_name_[minus/plus]_2.[fa/fq].\n");
    printf("  [--seed seed]: Optional argument, the seed used for simulation.\n");
    printf("  [--transcript name]: Optional argument, if set, only simulate reads from one transcript, the transcript name is 'name'.\n");
    printf("  [--read-model-file read_model_file]: Optional argument, use the read_model_file provided instead the one in stat folder of the 'sample_name'.\n");
    printf("  [--no-control]: Optional argument, indicate that there is no control.\n");
    return 0;
  }

  sprintf(refF, "%s.transcripts.fa", argv[1]);
  refs.readFrom(refF);
  M = refs.getM();
  
  N = atoi(argv[5]);
  seed = time(NULL);
  sim_tid = -1;
  readModelF[0] = 0;
  has_control = true;
  for (int i = 7; i < argc; ++i) {
    if (!strcmp(argv[i], "--seed")) {
      seed = 0;
      char *seedstr = argv[i + 1];
      int len = strlen(seedstr);
      for (int j = 0; j < len; ++j) seed = seed * 10 + (seedstr[j] - '0');
    }
    if (!strcmp(argv[i], "--transcript")) {
      string name = string(argv[i + 1]);
      for (int j = 1; j <= M; ++j) 
	if (refs.getRef(j)->getName() == name) { sim_tid = j; break; }
      if (sim_tid < 0) {
	fprintf(stderr, "Error: Cannot find transcript name %s from the reference!\n", argv[i + 1]);
	exit(-1);
      }
    }
    if (!strcmp(argv[i], "--read-model-file")) {
      strcpy(readModelF, argv[i + 1]);
    }
    if (!strcmp(argv[i], "--no-control")) has_control = false;
  }

  general_assert(has_control || !strcmp(argv[4], "plus"), "If no control data, only plus channel can be simulated!");

  // generate statName and readModelF (if necessary)
  string tmpStr = string(argv[3]);
  size_t strpos = tmpStr.find_last_of('/');
  if (strpos == string::npos) strpos = 0;
  else strpos = strpos + 1;
  assert(strpos < tmpStr.length());
  sprintf(statName, "%s.stat/%s", argv[3], tmpStr.substr(strpos).c_str());
  if (readModelF[0] == 0) {
    sprintf(readModelF, "%s_%s.read_model", statName, argv[4]);
  }


  sampler = new Sampler(seed);
  
  whole_model = new PROBerWholeModel(argv[2], (!strcmp(argv[4], "plus") && has_control ? 1 : 0), has_control);
  whole_model->read(argv[3], statName);

  read_model = new PROBerReadModel(&refs, sampler);
  read_model->read(readModelF);

  model_type = read_model->getModelType();

  if (model_type < 2) {
    sprintf(outF1, "%s_%s.%s", argv[6], argv[4], (model_type == 0 ? "fa" : "fq"));
    out1.open(outF1);
  }
  else {
    sprintf(outF1, "%s_%s_1.%s", argv[6], argv[4], (model_type == 2 ? "fa" : "fq"));
    sprintf(outF2, "%s_%s_2.%s", argv[6], argv[4], (model_type == 2 ? "fa" : "fq"));
    out1.open(outF1);
    out2.open(outF2);
  }

  int tid, pos, frag_len;
  whole_model->startSimulation(sim_tid);
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
