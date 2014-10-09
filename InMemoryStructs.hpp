#ifndef INMEMORYSTRUCTS_H_
#define INMEMORYSTRUCTS_H_

#include<vector>

// In memory alignment
struct InMemAlign {
  int tid, pos, fragment_length;
  double frac;
};

// In memory alignment group
struct InMemAlignG {
  int size; 
  InMemAlign **aligns; // alignments
  double noise_prob;

  ~InMemAlignG() {
    if (size > 0) {
      for (int i = 0; i < size; ++i) delete aligns[i];
      delete[] aligns;
    }
  }
};

#endif

