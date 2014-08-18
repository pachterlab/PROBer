#include<cstdio>
#include<cstring>
#include<cstdlib>
#include<cassert>
#include<string>
#include<vector>
#include<fstream>
#include<iostream>

#include "sam/bam.h"
#include "sam/sam.h"

#include "utils.h"
#include "my_assert.h"
#include "Transcripts.hpp"
#include "SEQstring.hpp"
#include "BamAlignment.hpp"
#include "AlignmentGroup.hpp"
#include "SamParser.hpp"
#include "BamWriter.hpp"
#include "MyHeap.hpp"

#include "RNASeqModel.hpp"

using namespace std;

int num_threads;
MyHeap my_heap; // a heap to record number of alignments contained in each partition

char tiF[STRLEN], bamOutF[STRLEN], cntF[STRLEN];
char paramsF[STRLEN], partitionF[STRLEN];

Transcripts transcripts;

const bam_header_t *header;

AlignmentGroup ag;
SamParser* parser;
vector<BamWriter*> writers;
BamWriter *writer0, *writer2;

READ_INT_TYPE N[4];
vector<READ_INT_TYPE> counts;
READ_INT_TYPE nUnique, nMulti, nIsoMulti;
HIT_INT_TYPE nHits;

bool bowtie_filter, polyA_filter;

RNASeqModel *model;

inline bool isGeneMultiRead(AlignmentGroup &ag) {
  int size = ag.size();
  if (size == 1) return false;
  string gene_id = transcripts.getTranscriptViaEid(ag.getAlignment(0)->getTid()).getGeneID();
  for (int i = 1; i < size; i++) 
    if (gene_id != transcripts.getTranscriptViaEid(ag.getAlignment(i)->getTid()).getGeneID()) return true;
  return false;
}

// In this version, only when all its alignments aligned to more than one isoform, the read is counted as an isoform multi-read 
inline bool isIsoMultiRead(AlignmentGroup &ag) {
  int size = ag.size();
  if (size == 1) return false;
  string iso_id = transcripts.getTranscriptViaEid(ag.getAlignment(0)->getTid()).getTranscriptID();
  for (int i = 1; i < size; i++) 
    if (iso_id != transcripts.getTranscriptViaEid(ag.getAlignment(i)->getTid()).getTranscriptID()) return true;
  return false;
}

// Filtering from bowtie
inline bool is_filtered_bowtie(AlignmentGroup &ag) {
  if (ag.isAligned()) return false;
  BamAlignment* ba = ag.getAlignment(0);

  uint8_t *p = NULL;
  char type = 0;

  if (ba->findTag("XM", p, type) && (type == 'i') && (ba->tag2i(p) > 0)) return true;
  if (ba->isPaired() && ba->findTag("XM", p, type, 2) && (type == 'i') && (ba->tag2i(p) > 0)) return true;
  return false;
}

// Filtering by sequence content
inline bool read_is_low_quality(AlignmentGroup &ag) {
  SEQstring si;
  int len, numA, numT;
  int threshold;

  ag.getSEQ(si);
  len = si.getLen();
  threshold = int(len * 0.8 + 0.5);
  numA = numT = 0;
  for (int i = 0; i < len; i++) { 
    if (si.baseAt(i) == 'A') ++numA;
    if (si.baseAt(i) == 'T') ++numT;
  }
  
  if ((numA >= threshold || numT >= threshold) && !ag.isPaired()) return true;
  
  ag.getSEQ(si, 2);
  len = si.getLen();
  threshold = int(len * 0.8 + 0.5);
  numA = numT = 0;
  for (int i = 0; i < len; i++) {
    if (si.baseAt(i) == 'A') ++numA;
    if (si.baseAt(i) == 'T') ++numT;
  }

  return (numA >= threshold || numT >= threshold);
}

inline bool filtered_due2polyA(AlignmentGroup &ag) {
  if (read_is_low_quality(ag)) return true;
  if (ag.isAligned()) {
    int s = ag.size();
    BamALignment *ba = NULL;
    for (int i = 0; i < s; ++i) {
      ba = ag.getAlignment(i);
      if (ba->getLeftMostPos() + MASK_LEN >= transcripts.getTranscriptViaEid(ba->getTid()).getLength()) {
	printf("Find a whatever alignment!\n");
	return true;
      }
    }     
  }
  return false;
}

void writeStat(const char* statName) {
  sprintf(cntF, "%s.cnt", statName);
  ofstream fout(cntF);
  assert(fout.is_open());
  N[3] = N[0] + N[1] + N[2];
  fout<< N[0]<< " "<< N[1]<< " "<< N[2]<< " "<< N[3]<< endl;
  fout<< nUnique<< " "<< nMulti<< " "<< nIsoMulti<< endl;
  fout<< nHits<< endl;

  // In this version, for all bins before the largest bin with at least one read, even if its bin count is 0, it is still printed out
  fout<< "0\t"<< N[0]<< endl;
  for (int i = 1; i < (int)counts.size(); i++) 
    fout<< i<< "\t"<< counts[i]<< endl;
  fout<< "Inf\t"<< N[2]<< endl;

  fout.close();
}

int main(int argc, char* argv[]) {
  if (argc < 7) { 
    printf("rsem-parse-alignments refName imdName statName number_of_partitions alignFType('s' for sam, 'b' for bam) alignF [-l fn_list] [--bowtie-filter] [--polyA-filter] [-q]\n");
    exit(-1);
  }

  // Load transcript information
  sprintf(tiF, "%s.ti", argv[1]);
  transcripts.readFrom(tiF);

  num_threads = atoi(argv[4]);
  assert(num_threads > 0);

  char *aux = NULL;
  bowtie_filter = polyA_filter = false;
  verbose = true;
  for (int i = 7; i < argc; i++) {
    if (!strcmp(argv[i], "-l")) aux = argv[i + 1];
    if (!strcmp(argv[i], "-q")) verbose = false;
    if (!strcmp(argv[i], "--bowtie-filter")) bowtie_filter = true;
    if (!strcmp(argv[i], "--polyA-filter")) polyA_filter = true;
  }

  parser = new SamParser(argv[5][0], argv[6], aux);
  header = parser->getHeader();
  transcripts.buildMappings(argv[2], header->n_targets, header->target_name);

  writers.assign(num_threads, NULL);
  writer0 = writer2 = NULL;
  for (int i = 0; i < num_threads; i++) {
    sprintf(bamOutF, "%s_%d.bam", argv[2], i);
    writers[i] = new BamWriter(bamOutF, (i == 0 ? header : NULL), "RSEM_intermediate");
  }
  sprintf(bamOutF, "%s_N0.bam", argv[2]);
  writer0 = new BamWriter(bamOutF, NULL, "RSEM_intermediate");
  sprintf(bamOutF, "%s_N2.bam", argv[2]);
  writer2 = new BamWriter(bamOutF, NULL, "RSEM_intermediate");

  sprintf(paramsF, "%s.mparams", argv[2]);
  model = new RNASeqModel(paramsF, 0);

  memset(N, 0, sizeof(N));
  counts.clear();
  nUnique = nMulti = nIsoMulti = 0;
  nHits = 0;

  READ_INT_TYPE cnt = 0;

  my_heap.init(num_threads);
  while (parser->next(ag)) {
    if ((polyA_filter && filtered_due2polyA(ag)) || (bowtie_filter && is_filtered_bowtie(ag))) {
      ++N[2];
      writer2->write(ag, 1);
    }
    else {
      bool isAligned = ag.isAligned();

      // update model parameters
      model->update_preprocess(ag, isAligned);

      if (isAligned) {
	++N[1];
	
	int id = my_heap.getTop();
	writers[id]->write(ag, 1); // remove seq and qual for secondary alignments
	my_heap.updateTop(ag.size());
	
	// Multi-read stats
	if (isGeneMultiRead(ag)) ++nMulti;
	else ++nUnique;
	if (isIsoMultiRead(ag)) ++nIsoMulti;
	
	nHits += (HIT_INT_TYPE)ag.size();
	if (ag.size() >= (int)counts.size()) counts.resize(ag.size() + 1, 0);
	++counts[ag.size()];
      }
      else {
	++N[0];
	writer0->write(ag, 1);
      }
      
      ++cnt;
      if (verbose && (cnt % 1000000 == 0)) cout<< cnt<< " reads are processed!"<< endl;
    }
  }

  sprintf(partitionF, "%s.partition", argv[2]);
  my_heap.print(partitionF);
  
  delete parser;
  for (int i = 0; i < num_threads; i++) delete writers[i];
  delete writer0;
  delete writer2;

  // release RNA-Seq Model
  model->write_preprocess(argv[2]);
  delete model;
  
  writeStat(argv[3]);

  return 0;
}
