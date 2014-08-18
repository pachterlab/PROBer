#include<cstdio>
#include<cassert>
#include<string>
#include<fstream>
#include<vector>

#include "utils.h"
#include "my_assert.h"
#include "RefSeqPolicy.h"
#include "PolyARules.h"
#include "RefSeq.hpp"
#include "Refs.hpp"

Refs::Refs() {
  M = 0;
  seqs.clear();
  has_polyA = false;
}

//inpF in fasta format
void Refs::makeRefs(char *inpF, RefSeqPolicy& policy, PolyARules& rules) {
  //read standard fasta format here
  std::ifstream fin;
  std::string tag, line, rawseq;

  seqs.clear();
  seqs.push_back(RefSeq()); // noise isoform

  M = 0;
  has_polyA = false;

  fin.open(inpF);
  general_assert(fin.is_open(), "Cannot open " + cstrtos(inpF) + "! It may not exist.");
  getline(fin, line);
  while ((fin) && (line[0] == '>')) {
    tag = line.substr(1);
    rawseq = "";
    while((getline(fin, line)) && (line[0] != '>')) {
      rawseq += line;
    }
    if (rawseq.size() <= 0) {
      printf("Warning: Fasta entry %s has an empty sequence! It is omitted!\n", tag.c_str()); 
      continue;
    }
    ++M;
    seqs.push_back(RefSeq(tag, policy.convert(rawseq), rules.getLenAt(tag)));
    has_polyA = has_polyA || seqs[M].getFullLen() < seqs[M].getTotLen();
  }
  fin.close();

  if (verbose) { printf("Refs.makeRefs finished!\n"); }
}

//inpF in fasta format, with sequence all in one line together
//option 0 read all, 1 do not read sequences
void Refs::loadRefs(char *inpF, int option) {
  std::ifstream fin;
  RefSeq seq;

  fin.open(inpF);
  general_assert(fin.is_open(), "Cannot open " + cstrtos(inpF) + "! It may not exist."); 
  seqs.clear();
  seqs.push_back(RefSeq());

  M = 0;
  has_polyA = false;

  bool success;
  do {
    success = seq.read(fin, option);
    if (success) {
    	seqs.push_back(seq);
        ++M;
    	has_polyA = has_polyA || seq.getFullLen() < seq.getTotLen();
    }
  } while (success);

  fin.close();

  assert(M + 1 == (int)seqs.size());

  if (verbose) { printf("Refs.loadRefs finished!\n"); }
}

void Refs::saveRefs(char* outF) {
  std::ofstream fout;

  fout.open(outF);
  for (int i = 1; i <= M; i++) {
    seqs[i].write(fout);
  }
  fout.close();

  if (verbose) { printf("Refs.saveRefs finished!\n"); }
}
