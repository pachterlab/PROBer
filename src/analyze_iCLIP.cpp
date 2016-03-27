#include<cstdio>
#include<cstring>
#include<cstdlib>
#include<cassert>
#include<string>
#include<map>
#include<vector>
#include<fstream>
#include<iostream>
#include<algorithm>
#include<pthread.h>
#include<unordered_map>

#include "htslib/sam.h"

#include "utils.h"
#include "my_assert.h"

#include "BamAlignment.hpp"
#include "AlignmentGroup.hpp"
#include "SamParser.hpp"
#include "BamWriter.hpp"

#include "PROBerReadModel_iCLIP.hpp"

#include "iCLIP_structs.hpp"
using namespace std;

bool verbose = true; // define verbose

// position, key for map
struct KeyType {
  char dir;
  int cid, pos;

  KeyType() : dir(0), cid(0), pos(0) {}
  
  KeyType(int cid, char dir, int pos) : dir(dir), cid(cid), pos(pos) {}
  
  bool operator< (const KeyType& o) const {
    if (cid != o.cid) return cid < o.cid;
    if (dir != o.dir) return dir < o.dir; // '+' < '-'
    return pos < o.pos;
  }
};

// read count weight, value for map
struct ValueType {
  double weight, ww; // weight, expected read counts at this position; ww, sum of weights in the window
  vector<double*> aligns; // point to expected weight from each alignment

  ValueType() : id(-1), weight(0), ww(0) { aligns.clear(); }
};

int n_pos; // number of distinct positions
map<KeyType, ValueType> posMap; // position map
map<KeyType, ValueType>::iterator iter; // iterator for posMap

int n_mhits; // number of multi-mapping reads' hits
double *fracs, *conprbs; // arrays storing multi-mapping reads expected weights and normalized sequencing error probabilities

// multi-mapping reads
struct MultiType {
  int offset, s, c; // offset, the offset at fracs and conprbs array; s, number of alignments; c, number of identical multi reads

  MultiType() : offset(0), s(0), c(0) {}

  MultiType(int offset, int s, int c = 1) : offset(offset), s(s), c(c) {}
};

int n_multi; // number of multi-mapping reads
vector<MultiType> multis; // multi-mapping reads

unordered_map<string, int> hash;


int model_type;

BamAlignment* ba;
AlignmentGroup ag;

PROBerReadModel_iCLIP* model;





/***  string variables  ***/

char imdName[STRLEN], statName[STRLEN];
char multiF[STRLEN], allF[STRLEN]; // multi-mapping reads, all reads
char modelF[STRLEN];

/***  optional arguments and auxiliary variables ***/  

bool bowtie_filter;
int max_hit_allowed; // maximum number of alignments allowed
int min_len; // minimum read length required
int max_len; // maximum read length
bool keep_alignments; // if keep the BAM file




/****************************************************************************************************/
// Parse alignments


// Filtering from bowtie
inline bool is_filtered_bowtie(AlignmentGroup &ag) {
  ba = ag.getAlignment(0);

  uint8_t *p = NULL;
  char type = 0;

  if (ba->findTag("XM", p, type) && (type == 'i') && (ba->tag2i(p) > 0)) return true;
  if (ba->isPaired() && ba->findTag("XM", p, type, 2) && (type == 'i') && (ba->tag2i(p) > 0)) return true;
  return false;
}

// Categorize reads and learn sequencing error model from uniquely-mapping reads
void parseAlignments(const char* alignF) {
  SamParser* parser;
  BamWriter *writer, *writerBam;

  READ_INT_TYPE N0, N11, N12, N2;

  pair<KeyType, ValueType> my_pair;
  
  parser = new SamParser(alignF);
  
  const char* program_id = parser->getProgramID();
  if (!strcmp(program_id, "Bowtie") || !strcmp(program_id, "bowtie")) bowtie_filter = true;
  
  const bam_hdr_t* header = parser->getHeader();

  sprintf(multiF, "%s_multi.bam", imdName);
  writer = new BamWriter(multiF, header, "PROBer iCLIP intermediate");
  if (keep_alignments) {
    sprintf(allF, "%s_alignments.bam", imdName);
    writerBam = new BamWriter(allF, header);
  }
    
  N0 = N11 = N12 = N2 = 0;
  n_mhits = 0;

  READ_INT_TYPE cnt = 0;
  
  while (parser->next(ag)) {
    if (keep_alignments) writerBam->write(ag);
    
    bool isAligned = ag.isAligned();

    if (ag.isFiltered() || ag.getSeqLength() < min_len || (ag.isPaired() && ag.getSeqLength(2) < min_len) || \
	(isAligned && ag.size() > max_hit_allowed) || (!isAligned && bowtie_filter && is_filtered_bowtie(ag))) {
      ++N2;
    }
    else if (isAligned) {
      // Read is alignable
      if (ag.size() == 1) {
	++N11;
	model->update(ag);

	ba = ag.getAlignment(0);
	my_pair.first.cid = ba->getTid();
	my_pair.first.dir = ba->getDir();
	my_pair.first.pos = ba->getPos();
	posMap.insert(my_pair).first->second.weight += 1.0;
      }
      else {
	++N12;
	n_mhits += ag.size();
	ag.sort_alignments();
	writer->write(ag);
      }
    }
    else {
      // Read is unalignable
      ++N0;
    }

    ++cnt;
    if (verbose && (cnt % 1000000 == 0)) cout<< cnt<< " reads are processed!"<< endl;
  }
  
  cout<< "N0 = "<< N0<< ", N11 = "<< N11<< ", N12 = "<< N12<< ", N2 = "<< N2<< ", n_mhits = "<< n_mhits<< endl;

  delete parser;
  delete writer;
  if (keep_alignments) delete writerBam;
}


/****************************************************************************************************/
// Estimate multi-mapping reads sequencing error probabilities

void processMultiReads() {
  SamParser *parser = new SamParser(multiF);

  pair<KeyType, ValueType> my_pair;

  
  ag.clear();

  fracs = new double[n_mhits];
  conprbs = new double[n_mhits];

  n_multi = 0; multis.clear();
  while (parser->next(ag)) {
  }

  delete parser;
  

  
  n_multi = multis.size();
  for (int i = 0; i < n_multi; ++i) {
    assert(parser->next(ag));
    model->calcProbs(ag, multis[i].aligns);
    if (verbose && (i + 1) % 1000000 == 0) cout<< i + 1<< " multi-mapping reads are processed!"<< endl;
  }
  delete parser;

  cout<< "Setting up probabilities is done!"<< endl;

  for (int i = 0; i < n_multi; ++i) multis[i].sort_alignments();
  sort(multis.begin(), multis.end());

  cout<< "Sorting is finished."<< endl;

  char tmpF[STRLEN];
  sprintf(tmpF, "%s.tmp", imdName);
  FILE *fo = fopen(tmpF, "w");
  
  int pos = 0;
  for (int i = 1; i < n_multi; ++i)
    if (multis[pos] != multis[i]) {
      if (i - pos > 1) {
	for (int j = pos; j < i; ++j) {
	  for (int k = 0; k < multis[j].s; ++k) fprintf(fo, "%.6g ", multis[j].aligns[k].conprb);
	  fprintf(fo, "%c", j == i - 1 ? '\n' : '\t');
	}
      }
      pos = i;
    }
  if (n_multi - pos > 1) {
    for (int j = pos; j < n_multi; ++j) {
      for (int k = 0; k < multis[j].s; ++k) fprintf(fo, "%.6g ", multis[j].aligns[k].conprb);
      fprintf(fo, "%c", j == n_multi - 1 ? '\n' : '\t');
    }
  }
  fclose(fo);
  cout<< "Done!"<< endl;
}


/****************************************************************************************************/


void init() {
  model = new PROBerReadModel_iCLIP(model_type, max_len);

  n_pos = 0; posMap.clear();

  fracs = conprbs = NULL;
}

void release() {
  sprintf(modelF, "%s.model", imdName);
  model->write(modelF);
  delete model;

  delete[] fracs;
  delete[] conprbs;
}


/****************************************************************************************************/


int main(int argc, char* argv[]) {
  if (argc < 4) { 
    printf("PROBer-analyze-iCLIP model_type imdName alignF [-m max_hit_allowed][--shorter-than min_len] [--keep-alignments] [--max-len max_len] [-q]\n");
    exit(-1);
  }

  model_type = atoi(argv[1]);
  strcpy(imdName, argv[2]);
  
  bowtie_filter = false;
  max_hit_allowed = 2147483647; // 2^31 - 1
  min_len = -1;
  max_len = -1;
  keep_alignments = false;
  
  for (int i = 4; i < argc; i++) {
    if (!strcmp(argv[i], "-q")) verbose = false;
    if (!strcmp(argv[i], "-m")) max_hit_allowed = atoi(argv[i + 1]);
    if (!strcmp(argv[i], "--shorter-than")) min_len = atoi(argv[i + 1]);
    if (!strcmp(argv[i], "--keep-alignments")) keep_alignments = true;
    if (!strcmp(argv[i], "--max-len")) max_len = atoi(argv[i + 1]);
  }

  init();  
  parseAlignments(argv[3]);
  model->finish();
  processMultiReads();
  release();

  return 0;
}
