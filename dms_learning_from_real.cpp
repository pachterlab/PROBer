#include<cstdio>
#include<cstring>
#include<cstdlib>
#include<cassert>

#include "SamParser.hpp"
#include "BamAlignment.hpp"
#include "AlignmentGroup.hpp"
#include "DMSWholeModel.hpp"

using namespace std;

int num_threads;

SamParser *in;
AlignmentGroup ag;

DMSWholeModel *model;

int choice;
double value;

int main(int argc, char* argv[]) {
  if (argc < 6 || argc > 8) {
    printf("Usage: dms_learning_from_real config_file input.bam output_name num_threads <eXpress (0) or RSEM (1)> [--gamma input_file_name]\n");
    exit(-1);
  }

  num_threads = atoi(argv[4]);
  choice = atoi(argv[5]);

  in = new SamParser('b', argv[2]);
  model = new DMSWholeModel(argv[1], in->getHeader(), num_threads);

  int cnt = 0;
  int wrong_dir = 0;
  // reading counts
  model->init();
  while (in->next(ag)) {
    ++cnt;
    if (!ag.isAligned()) model->update(0, 0, 1.0);
    else {
      int s = ag.size();
      double left = 1.0;
      for (int i = 0; i < s; i++) {
	BamAlignment *ba = ag.getAlignment(i);
	if (ba->getDir() == '+') { 
	  if (choice == 0) { 
	    uint8_t *p;
	    char type;
	    assert(ba->findTag("XP", p, type) && type == 'f');
	    value = ba->tag2f(p);
	  }
	  else {
	    value = ba->getFrac();
	    assert(value >= 0.0);
	  }
	  model->update(ba->getTid() + 1, ba->getLeftMostPos(), value);
	  left -= value;
	}
	else ++wrong_dir;
      }

      if (left > 1e-8) {
	assert(choice == 1);
	model->update(0, 0, left);
      }
    }      
    if (cnt % 1000000 == 0) printf("%d reads loaded!\n", cnt);
  }

  printf("WRONG_DIR = %d\n", wrong_dir);

  if (argc == 8) {
    assert(!strcmp(argv[6], "--gamma"));
    model->read(argv[7]);
  }

  model->runEM(200);
  model->write(argv[3]);

  delete in;
  delete model;

  return 0;
}
