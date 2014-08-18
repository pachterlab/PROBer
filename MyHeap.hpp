#ifndef MYHEAP_H_
#define MYHEAP_H_

#include<cassert>
#include<vector>
#include<fstream>
#include<iostream>

#include "utils.h"

struct HeapType {
  int id;
  READ_INT_TYPE nreads;
  HIT_INT_TYPE nlines;
  
  HeapType() { id = -1; nreads = 0; nlines = 0; }
  HeapType(int id, READ_INT_TYPE nreads, HIT_INT_TYPE nlines) { this->id = id; this->nreads = nreads; this->nlines = nlines; }

  HeapType(const HeapType& o) { id = o.id; nreads = o.nreads; nlines = o.nlines; }

  bool operator< (const HeapType& o) const {
    return nlines < o.nlines || (nlines == o.nlines && id < o.id);
  }
};

class MyHeap {
public:
  MyHeap() { size = 0; elements.clear(); }

  void init(int size);
  int getTop() const { return elements[0].id; }
  void updateTop(HIT_INT_TYPE nlines);
  
  void print(const char* outF = NULL);

private:
  int size;
  std::vector<HeapType> elements;
};

inline void MyHeap::init(int size) {
  this->size = size;
  elements.clear();
  for (int i = 0; i < size; i++) 
    elements.push_back(HeapType(i, 0, 0));
}

inline void MyHeap::updateTop(HIT_INT_TYPE nlines) {
  int pos = 0, minpos;
  ++elements[pos].nreads;
  elements[pos].nlines += nlines;
  while (pos * 2 + 1 < size) {
    minpos = pos; 
    if (elements[pos * 2 + 1] < elements[pos]) minpos = pos * 2 + 1;
    if ((pos * 2 + 2 < size) && (elements[pos * 2 + 2] < elements[minpos])) minpos = pos * 2 + 2;
    if (pos == minpos) break;
    std::swap(elements[pos], elements[minpos]);
    pos = minpos;
  }
}

inline void MyHeap::print(const char* outF) {
  if (outF != NULL) {
    std::ofstream fout(outF);
    assert(fout.is_open());
    
    for (int i = 0; i < size; i++) 
      fout<< elements[i].id<< '\t'<< elements[i].nreads<< '\t'<< elements[i].nlines<< std::endl;
  }

  if (verbose) {
    for (int i = 0; i < size; i++) 
      std::cout<< "Partition "<< elements[i].id<< " has "<< elements[i].nreads<< " reads and "<< elements[i].nlines<< " lines."<< std::endl;
  }
}

#endif
