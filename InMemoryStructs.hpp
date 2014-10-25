#ifndef INMEMORYSTRUCTS_H_
#define INMEMORYSTRUCTS_H_

#include<vector>

#include "utils.h"

// In memory alignment
struct InMemAlign {
  int tid, pos, fragment_length;
  double frac; // the conditional probability of generating the read based on read_model; After the main EM, the expected weight will be assigned to this field

  InMemAlign() : tid(0), pos(0), fragment_length(0), frac(0.0) {}
};

// In memory alignment group
struct InMemAlignG {
  int size; 
  double noise_prob;

  static InMemAlign* aligns; // alignments

  InMemAlignG() : size(0), noise_prob(0.0) {}
};

// Store in memory information for all alignments of a thread
struct InMemChunk {
  READ_INT_TYPE pos, nreads;
  InMemAlign *aligns;
  InMemAlignG *reads;

  InMemChunk(READ_INT_TYPE nreads, HIT_INT_TYPE nlines) {
    this->nreads = nreads;
    reads = new InMemAlignG[nreads];
    --reads; // because reads[1] is the start position
    aligns = new InMemAlign[nlines];

    reset();
  }

  /*
    @func   reset InMemAlignG.aligns to NULL
   */
  void reset() {
    pos = 0;
    InMemAlignG::aligns = aligns;
  }

  /*
    @func  return the read current pointer points to
   */
  InMemAlignG* next() {
    if (pos >= nreads) return NULL;
    if (pos > 0) InMemAlignG::aligns += reads[pos].size;
    ++pos;
    return reads + pos;
  }

  ~InMemChunk() {
    ++reads;
    delete[] reads;
    delete[] aligns;
  }  
};

#endif

