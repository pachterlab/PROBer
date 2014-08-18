#ifndef TRANSALIGNMENTGROUP_H_
#define TRANSALIGNMENTGROUP_H_

#include<cstdio>
#include<cstring>
#include<cstdlib>
#include<cassert>
#include<map>
#include<vector>
#include<string>
#include<algorithm>

#include "sam/sam.h"

#include "Transcript.hpp"
#include "Transcripts.hpp"
#include "TransBamAlignment.hpp"
#include "TransAlignmentGroup.hpp"

class TransAlignmentGroup {
public:
  TransAlignmentGroup(const Transcripts& transcripts, const std::map<std::string, int>& refmap) : transcripts(transcripts), refmap(refmap) { 
    leftover = -1;
    size = max_size = 0;
    alignments.clear(); 
  }
  
  TransAlignmentGroup(const TransAlignmentGroup& o) : transcripts(o.transcripts), refmap(o.refmap), leftover(o.leftover), size(o.size), max_size(o.max_size), alignments(o.alignments) {
    assert(o.max_size == 0); // temporarily put this assert here 
  }

  ~TransAlignmentGroup() {
    for (int i = 0; i < max_size; i++) delete alignments[i];
  }

  void allocate() {
    if (size >= max_size) { 
      alignments.push_back(new TransBamAlignment());
      ++max_size;
    }    
  }

  bool read(samfile_t*);
  bool convert(samfile_t*); // convert and then write down the results

private:
  const Transcripts& transcripts;
  const std::map<std::string, int>& refmap;

  char leftover;
  int size, max_size;
  std::vector<TransBamAlignment*> alignments;

  struct doCompare {
    const TransAlignmentGroup& info;
    
    doCompare(const TransAlignmentGroup& info) : info(info) {}
    
    // test if lhs < rhs
    bool operator() (const int& lhs, const int& rhs) const {
      return info.alignments[lhs]->compare(info.alignments[rhs]) < 0;
    }
  };
};

inline bool TransAlignmentGroup::read(samfile_t *in) {
  const char* cname = NULL;
  TransBamAlignment *tmp = NULL;

  switch (leftover) {
  case 1 : 
    // swap
    tmp = alignments[size]; alignments[size] = alignments[0]; alignments[0] = tmp;
    break;
  case 0 : return false;
  case -1 : 
    size = 0; allocate(); 
    leftover = alignments[size]->read(in);
    if (leftover == 0) return false;
    break;
  default : assert(false);
  }
  
  cname = alignments[0]->getName();
  assert(cname[0] != 0);
  size = 1;

  while ((allocate(), (leftover = alignments[size]->read(in))) && (!strcmp(cname, alignments[size]->getName()))) {
    assert(alignments[size]->isPaired() == alignments[0]->isPaired());
    ++size;
  }

  return true;
}

inline bool TransAlignmentGroup::convert(samfile_t *out) {
  assert(size > 0);
  if (alignments[0]->isAligned() == 0) {
    assert(size == 1);
    alignments[0]->write(out);
    return true;
  }

  // Convert each transcript alignment
  std::map<std::string, int>::const_iterator iter;
  for (int i = 0; i < size; i++) {
    const Transcript& transcript = transcripts.getTranscriptViaEid(alignments[i]->getTid());
    iter = refmap.find(transcript.getSeqName());
    assert(iter != refmap.end());
    alignments[i]->convert(iter->second, transcript);
  }

  if (size == 1) { 
    alignments[0]->write(out);
    return true;
  }

  // Sort and merge alignments
  std::vector<int> ids(size, -1);
  for (int i = 0; i < size; i++) ids[i] = i;
  sort(ids.begin(), ids.end(), doCompare(*this));
  
  int curid = 0;
  double frac = alignments[ids[curid]]->getFrac();
  for (int i = 1; i < size; i++) {
    if (alignments[ids[curid]]->compare(alignments[ids[i]]) != 0) {
      if (frac >= 0.0) alignments[ids[curid]]->setFrac(frac);
      alignments[ids[curid]]->write(out);
      curid = i; frac = 0.0;
    }
    if (frac < 0.0) assert(alignments[ids[i]]->getFrac() < 0.0);
    else frac += alignments[ids[i]]->getFrac();
  }
  if (frac >= 0.0) alignments[ids[curid]]->setFrac(frac);
  alignments[ids[curid]]->write(out);

  return true;
}

#endif
