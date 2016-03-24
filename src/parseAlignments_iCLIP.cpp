#include<cstdio>
#include<cstring>
#include<cstdlib>
#include<cassert>
#include<string>
#include<vector>
#include<fstream>
#include<iostream>

#include "htslib/sam.h"

#include "utils.h"
#include "my_assert.h"
#include "BamAlignment.hpp"
#include "AlignmentGroup.hpp"
#include "SamParser.hpp"
#include "BamWriter.hpp"

using namespace std;

bool verbose = true; // define verbose

char imdName[STRLEN], statName[STRLEN];
char bgF[STRLEN], uniqF[STRLEN], multiF[STRLEN], filtF[STRLEN]; // background, unique, multi-mapping, filtered

const bam_hdr_t *header;

AlignmentGroup ag;
SamParser* parser;
BamWriter *writer0, *writer11, *writer12, *writer2;

READ_INT_TYPE N0, N11, N12, N2;
HIT_INT_TYPE nHits;

bool bowtie_filter;
int max_hit_allowed; // maximum number of alignments allowed
int min_len; // minimum read length required

// Filtering from bowtie
inline bool is_filtered_bowtie(AlignmentGroup &ag) {
  BamAlignment* ba = ag.getAlignment(0);

  uint8_t *p = NULL;
  char type = 0;

  if (ba->findTag("XM", p, type) && (type == 'i') && (ba->tag2i(p) > 0)) return true;
  if (ba->isPaired() && ba->findTag("XM", p, type, 2) && (type == 'i') && (ba->tag2i(p) > 0)) return true;
  return false;
}

int main(int argc, char* argv[]) {
  if (argc < 3) { 
    printf("PROBer-parse-alignments-iCLIP imdName alignF [-m max_hit_allowed][--shorter-than min_len] [-q]\n");
    exit(-1);
  }

  bowtie_filter = false;
  max_hit_allowed = 2147483647; // 2^31 - 1
  min_len = -1;

  for (int i = 3; i < argc; i++) {
    if (!strcmp(argv[i], "-q")) verbose = false;
    if (!strcmp(argv[i], "-m")) max_hit_allowed = atoi(argv[i + 1]);
    if (!strcmp(argv[i], "--shorter-than")) min_len = atoi(argv[i + 1]);
  }

  parser = new SamParser(argv[2]);

  const char* program_id = parser->getProgramID();
  if (!strcmp(program_id, "Bowtie") || !strcmp(program_id, "bowtie")) bowtie_filter = true;

  header = parser->getHeader();

  sprintf(bgF, "%s_bg.bam", imdName);
  writer0 = new BamWriter(bgF, header, "PROBer iCLIP intermediate");
  sprintf(uniqF, "%s_uniq.bam", imdName);
  writer11 = new BamWriter(uniqF, header, "PROBer iCLIP intermediate");
  sprintf(multiF, "%s_multi.bam", imdName);
  writer12 = new BamWriter(multiF, header, "PROBer iCLIP intermediate");
  sprintf(filtF, "%s_filt.bam", imdName);
  writer2 = new BamWriter(filtF, header, "PROBer iCLIP intermediate");

  N0 = N11 = N12 = N2 = 0;
  nHits = 0;

  READ_INT_TYPE cnt = 0;

  while (parser->next(ag)) {
    bool isAligned = ag.isAligned();

    if (ag.isFiltered() || ag.getSeqLength() < min_len || (ag.isPaired() && ag.getSeqLength(2) < min_len) || \
	(isAligned && ag.size() > max_hit_allowed) || (!isAligned && bowtie_filter && is_filtered_bowtie(ag))) {
      ++N2;
      writer2->write(ag);
    }
    else if (isAligned) {
      // Read is alignable
      if (ag.size() == 1) {
	++N11;
	writer11->write(ag);
      }
      else {
	++N12;
	writer12->write(ag);
      }
      
      nHits += (HIT_INT_TYPE)ag.size();
    }
    else {
      // Read is unalignable
      ++N0;
      writer0->write(ag);
    }

    ++cnt;
    if (verbose && (cnt % 1000000 == 0)) cout<< cnt<< " reads are processed!"<< endl;
  }

  printf("N0 = %d, N11 = %d, N12 = %d, N2 = %d\n", N0, N11, N12, N2);
  
  delete parser;
  delete writer0;
  delete writer11;
  delete writer12;
  delete writer2;

  return 0;
}
