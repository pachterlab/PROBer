#include<cstdio>
#include<cstring>
#include<cstdlib>
#include<cassert>

#include "sam/bam.h"
#include "DMSWholeModel.hpp"

using namespace std;

int num_threads;
bam_header_t *header;
DMSWholeModel *model;

int main(int argc, char* argv[]) {
  if (argc < 7 || argc > 9) {
    printf("Usage: dms_learning_from_simulated config_file header_file input_sim_file output_name num_threads <'SE' or 'PE'> [--gamma input_file_name]\n");
    exit(-1);
  }

  num_threads = atoi(argv[5]);
  header = sam_header_read2(argv[2]);
  model = new DMSWholeModel(argv[1], header, num_threads);

  if (argc == 9) {
    assert(!strcmp(argv[7], "--gamma"));
    model->read(argv[8]);
  } 

  FILE *fi = fopen(argv[3], "r");
  model->init();
  int tid, pos, frag_len;
  int cnt = 0;
  while (fscanf(fi, "%d %d %d", &tid, &pos, &frag_len) == 3) {
    if (!strcmp(argv[6], "SE")) model->update(tid, pos, 1.0);
    else model->update(tid, pos, frag_len, 1.0);
    ++cnt;
    if (cnt % 1000000 == 0) printf("%d updated!\n", cnt);
  }
  fclose(fi);

  model->runEM(200);
  model->write(argv[4]);

  delete header;
  delete model;

  return 0;
}
