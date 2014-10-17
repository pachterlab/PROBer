#ifndef INMEMORYSTRUCTS_H_
#define INMEMORYSTRUCTS_H_

#include<vector>

// In memory alignment
struct InMemAlign {
  int tid, pos, fragment_length;
  double frac; // the conditional probability of generating the read based on read_model

  InMemAlign(int tid, int pos, int fragment_length, double frac) : tid(tid), pos(pos), fragment_length(fragment_length), frac(frac) {
  }
};

// In memory alignment group
struct InMemAlignG {
  int size; 
  InMemAlign **aligns; // alignments
  double noise_prob;

  InMemAlignG() {
    size = 0; 
    aligns = NULL;
    noise_prob = 0.0;
  }

  ~InMemAlignG() {
    if (size > 0) {
      for (int i = 0; i < size; ++i) delete aligns[i];
      delete[] aligns;
    }
  }
};

#endif

