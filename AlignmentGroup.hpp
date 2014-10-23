#ifndef ALIGNMENTGROUP_H_
#define ALIGNMENTGROUP_H_

#include<cstdio>
#include<cstring>
#include<cstdlib>
#include<cassert>
#include<vector>

#include "sam/sam.h"

#include "SEQstring.hpp"
#include "QUALstring.hpp"
#include "BamAlignment.hpp"

class AlignmentGroup {
public:
  AlignmentGroup() { 
    leftover = -1;
    s = max_size = 0;
    alignments.clear(); 
  }

  AlignmentGroup(const AlignmentGroup& o) : leftover(o.leftover), s(o.s), max_size(o.max_size), alignments(o.alignments) {
    assert(o.max_size == 0); // temporarily put this assert here
  }

  ~AlignmentGroup() {
    for (int i = 0; i < max_size; i++) delete alignments[i];
  }

  void allocate() {
    if (s >= max_size) { 
      alignments.push_back(new BamAlignment());
      ++max_size;
    }    
  }

  void markAsFiltered() {
    for (int i = 0; i < s; ++i) alignments[i]->markAsFiltered();
  }

  bool read(samfile_t*);
  bool write(samfile_t*, int=0);

  bool isPaired() const { return (s > 0) && alignments[0]->isPaired(); }

  bool isAligned() const { return (s > 0) && (alignments[0]->isAligned() > 0); }

  int size() const { return s; }

  std::string getName(int mate = 0) const { 
    assert(s > 0);
    return alignments[0]->getName(mate);
  }
  
  int getSeqLength(int mate = 1) const { 
    assert(s > 0);
    return alignments[0]->getSeqLength(mate);
  }

  bool getSEQ(SEQstring& si, int mate = 1) {
    assert(s > 0); 
    return alignments[0]->getSEQ(si, mate);
  }

  bool getQUAL(QUALstring& qi, int mate = 1) {
    assert(s > 0);
    return alignments[0]->getQUAL(qi, mate);
  }

  BamAlignment* getAlignment(int id) { 
    assert(id >=0 && id < s);
    return alignments[id];
  }

private:
  char leftover;
  int s, max_size;
  std::vector<BamAlignment*> alignments;
};

inline bool AlignmentGroup::read(samfile_t *in) {
  const char *cname = NULL, *name = NULL;
  BamAlignment *tmp = NULL;

  switch (leftover) {
  case 1 :
    // swap
    tmp = alignments[s]; alignments[s] = alignments[0]; alignments[0] = tmp;
    break;
  case 0 : return false;
  case -1 :
    s = 0; allocate();
    leftover = alignments[s]->read(in);
    if (leftover == 0) return false;
    break;
  default : assert(false);
  }

  cname = alignments[0]->getName();
  assert(cname[0] != 0);
  s = 1;

  while (allocate(), (leftover = alignments[s]->read(in, alignments[0]))) {
    name = alignments[s]->getName();
    if (name[0] != 0 && strcmp(cname, name)) break;
    assert(alignments[s]->isPaired() == alignments[0]->isPaired());
    ++s;
  }

  return true;
}

// choice: 0, do nothing; 1, delete read sequence and qual score; 2, add read sequence and qual score.
inline bool AlignmentGroup::write(samfile_t *out, int choice) {
  assert(s > 0);
  alignments[0]->write(out);
  for (int i = 1; i < s; i++) alignments[i]->write(out, choice, alignments[0]);
  return true;
}

#endif
