#ifndef INMEMORYSTRUCTS_H_
#define INMEMORYSTRUCTS_H_

#include<vector>

#include "utils.h"

// In memory alignment
struct InMemAlign {
  int tid, pos, fragment_length;
  double conprb, frac; // conprb, the conditional probability of generating the read based on read_model; frac the expected weight

  InMemAlign() : tid(0), pos(0), fragment_length(0), conprb(0.0), frac(0.0) {}
};

// In memory alignment group
struct InMemAlignG {
  int size; 
  double noise_conprb;

  InMemAlignG() : size(0), noise_conprb(0.0) {}
};

// Store in memory information for all alignments of a thread
struct InMemChunk {
  READ_INT_TYPE pos, nreads;
  InMemAlign *aligns, *pointer;
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
    pointer = aligns;
  }

  /*
    @param   aread     In memory read group information
    @param   alignArr  In memory alignment pointer
    @return  true if has more reads, false otherwise
   */
  bool next(InMemAlignG*& aread, InMemAlign*& alignArr) {
    if (pos >= nreads) return false;
    if (pos > 0) pointer += reads[pos].size;
    ++pos;

    aread = reads + pos;
    alignArr = pointer;

    return true;
  }

  ~InMemChunk() {
    ++reads;
    delete[] reads;
    delete[] aligns;
  }  
};

#endif

