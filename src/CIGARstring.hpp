#ifndef CIGARSTRING_H_
#define CIGARSTRING_H_

#include<cassert>
#include<string>
#include<sstream>

#include<stdint.h>
#include "sam/bam.h"

// An iterator for BAM CIGAR strings
class CIGARstring {
public:
  CIGARstring() : cigar(NULL), dir(0), len(0) {}
  
  void setUp(const uint32_t *cigar, int len, char dir = '+') { 
    this->cigar = cigar;
    this->len = len;
    
    assert(dir == '+' || dir == '-');
    this->dir = dir;
  }

  void setDir(char dir) { 
    assert(dir == '+' || dir == '-');
    this->dir = dir;
  }

  int getLen() const { return len; }

  int opAt(int pos) const {
    assert(pos >= 0 && pos < len);
    return bam_cigar_op(dir == '+' ? cigar[pos] : cigar[len - pos - 1]);
  }

  char opchrAt(int pos) const { 
    assert(pos >= 0 && pos < len);
    return bam_cigar_opchr(dir == '+' ? cigar[pos] : cigar[len - pos - 1]);
  }

  int oplenAt(int pos) const {
    assert(pos >= 0 && pos < len);
    return bam_cigar_oplen(dir == '+' ? cigar[pos] : cigar[len - pos - 1]);
  }

  // toString will reset dir
  std::string toString(char dir = '+') {
    assert(dir == '+' || dir == '-');
    this->dir = dir;
    std::ostringstream strout;
    for (int i = 0; i < len; i++) strout<< opAt(i)<< oplenAt(i);
    this->dir = 0;
    return strout.str();
  }

private:
  const uint32_t *cigar;
  char dir;
  int len; // len, cigar string length
};

#endif
