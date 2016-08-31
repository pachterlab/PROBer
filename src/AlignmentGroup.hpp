#ifndef ALIGNMENTGROUP_H_
#define ALIGNMENTGROUP_H_

#include<cstdio>
#include<cstring>
#include<cstdlib>
#include<cassert>
#include<vector>
#include<algorithm>

#include "htslib/sam.h"

#include "SEQstring.hpp"
#include "QUALstring.hpp"
#include "BamAlignment.hpp"


// comparison function for iCLIP
inline bool cmp_iCLIP(const BamAlignment* a, const BamAlignment* b) {
  if (a->getTid() != b->getTid()) return a->getTid() < b->getTid();
  if (a->getDir() != b->getDir()) return a->getDir() < b->getDir();
  return a->getCrosslinkSite() < b->getCrosslinkSite();
}

// One alignment group associates with one file. 
class AlignmentGroup {
public:
  AlignmentGroup() { 
    leftover = -1;
    s = max_size = 0;
    alignments.clear(); 
  }

  ~AlignmentGroup() {
    for (int i = 0; i < max_size; ++i) delete alignments[i];
  }

  /*
    @func   clear the alignment group so that we can use it for next SAM/BAM file
   */
  void clear() {
    leftover = -1;
    s = 0;
  }

  bool isFiltered() const { return (s > 0) && alignments[0]->isFiltered(); }

  void markAsFiltered() {
    for (int i = 0; i < s; ++i) alignments[i]->markAsFiltered();
  }

  bool read(samFile* in, bam_hdr_t* header);
  bool write(samFile* out, bam_hdr_t* header, int choice = 0); // only writ out one alignment for filtered reads

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

  void sort_alignments() {
    assert(s > 1);
    std::sort(alignments.begin(), alignments.begin() + s, cmp_iCLIP);
  }
  
private:
  char leftover; // if has next read. -1, initial value, stands for not called; 0, no next read; 1, has next read, which is alignments[s]
  int s, max_size; // s, total number of alignments; max_size, max capacity
  std::vector<BamAlignment*> alignments; // pointers to BamAlignment objects

  void allocate() {
    if (s >= max_size) { 
      alignments.push_back(new BamAlignment());
      ++max_size;
    }    
  }
};

inline bool AlignmentGroup::read(samFile* in, bam_hdr_t* header) {
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
    leftover = alignments[s]->read(in, header);
    if (leftover == 0) return false;
    break;
  default : assert(false);
  }

  cname = alignments[0]->getName();
  assert(cname[0] != 0);
  s = 1;

  while (allocate(), (leftover = alignments[s]->read(in, header, alignments[0]))) {
    name = alignments[s]->getName();
    if (name[0] != 0 && strcmp(cname, name)) break;
    assert(alignments[s]->isPaired() == alignments[0]->isPaired());
    ++s;
  }

  return true;
}

// choice: 0, do nothing; 1, delete read sequence and qual score; 2, add read sequence and qual score.
inline bool AlignmentGroup::write(samFile *out, bam_hdr_t* header, int choice) {
  assert(s > 0);
  alignments[0]->write(out, header);
  if (s > 1 && !alignments[0]->isFiltered())
    for (int i = 1; i < s; ++i) alignments[i]->write(out, header, choice, alignments[0]);
  return true;
}

#endif
