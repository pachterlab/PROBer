#include<ctime>
#include<cstdio>
#include<cstring>
#include<cstdlib>
#include<cassert>

#include "sampling.hpp"
#include "DMSWholeModel.hpp"

using namespace std;

seedType seed;
Sampler *sampler;
DMSWholeModel *model;

int num_reads;
FILE *fo;
char outF[1005];

int main(int argc, char* argv[]) {
  if (argc < 6 || argc > 7) {
    printf("Usage: dms_simulate_reads config_file input_name output_name <'minus' or 'plus'> number_of_reads [seed]\n");
    exit(-1);
  }

  num_reads = atoi(argv[5]);
  if (argc == 7) {
    int len = strlen(argv[6]);
    seed = 0;
    for (int i = 0; i < len; ++i) seed = seed * 10 + (argv[6][i] - '0');
  }
  else seed = time(NULL);

  sampler = new Sampler(seed);

  model = new DMSWholeModel(argv[1]);
  model->read(argv[2]);
  if (!strcmp(argv[4], "plus")) model->read(argv[2]);

  sprintf(outF, "%s_%s.txt", argv[3], argv[4]);
  fo = fopen(outF, "w");

  int tid, pos, frag_len;

  model->startSimulation();
  for (int i = 0; i < num_reads; ++i) {
    model->simulate(sampler, tid, pos, frag_len);
    fprintf(fo, "%d\t%d\t%d\n", tid, pos, frag_len);
    if ((i + 1) % 1000000 == 0) printf("%d FIN!\n", i + 1);
  }
  model->finishSimulation();

  fclose(fo);

  delete model;

  return 0;
}
