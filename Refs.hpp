#ifndef REFS_H_
#define REFS_H_

#include<vector>

#include "RefSeqPolicy.h"
#include "PolyARules.h"
#include "RefSeq.hpp"

class Refs {
 public:
  Refs();

  void makeRefs(char*, RefSeqPolicy&, PolyARules&);
  void loadRefs(char*, int = 0);
  void saveRefs(char*);

  int getM() { return M; } // get number of isoforms

  RefSeq& getRef(int sid) { return seqs[sid]; } // get a particular reference

  bool hasPolyA() { return has_polyA; } // if any of sequence has poly(A) tail added

 private:
  int M; // # of isoforms, id starts from 1
  std::vector<RefSeq> seqs;  // reference sequences, starts from 1; 0 is for noise gene
  bool has_polyA; // if at least one sequence has polyA added, the value is true; otherwise, the value is false
};

#endif
