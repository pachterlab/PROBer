#include<cstdio>
#include<cstring>
#include<cstdlib>
#include<cassert>

#include "stdint.h"
#include "sam/bam.h"
#include "sam/sam.h"

using namespace std;

samfile_t *in;
bam_header_t *header;
bam1_t *b;

int tid;
int len;

double *counts;

int main(int argc, char* argv[]) {
  if (argc != 4) {
    printf("Usage: extract_count_vector isoform_name input.bam output.txt\n");
    exit(-1);
  }

  in = samopen(argv[2], "rb", NULL);
  header = in->header;
  
  tid = -1;
  for (int i = 0; i < header->n_targets; ++i) {
    if (!strcmp(header->target_name[i], argv[1])) {
      tid = i;
      len = header->target_len[i];
      break;
    }
  }
  assert(tid >= 0);

  counts = new double[len];
  memset(counts, 0, sizeof(double) * len);

  int cnt = 0;

  b = bam_init1();
  while (samread(in, b) > 0) {
    if (!(b->core.flag & 0x0004) && b->core.tid == tid) {
      assert(!(b->core.flag & 0x0010));
      uint8_t *p_tag = bam_aux_get(b, "XP");
      assert(p_tag != NULL);
      counts[b->core.pos] += bam_aux2f(p_tag);
    }
    ++cnt;
    if (cnt % 1000000 == 0) printf("FIN %d\n", cnt);
  }

  FILE *fo = fopen(argv[3], "w");
  fprintf(fo, "%d", len);
  for (int i = 0; i < len; ++i) fprintf(fo, "\t%.2f", counts[i]);
  fprintf(fo, "\n");
  fclose(fo);

  delete[] counts;

  return 0;
}
