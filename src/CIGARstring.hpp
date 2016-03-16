#ifndef CIGARSTRING_H_
#define CIGARSTRING_H_

#include<cassert>
#include<string>
#include<sstream>

#include<stdint.h>
#include "htslib/sam.h"

// An iterator for BAM CIGAR strings
class CIGARstring {
public:
  CIGARstring() : cigar(NULL), is_ori(true), return_current(true), len(0) {}
  
  void setUp(const uint32_t *cigar, int len, bool is_ori) { 
    this->cigar = cigar;
    this->len = len;
    this->is_ori = is_ori;
    // Default, return the CIGAR string in original read sequence's order
    return_current = (is_ori ? true : false);
  }
  
  // '+' returns CIGAR string in original read sequence's order; '-' returns the reverse CIGAR
  void setDir(char dir) { 
    assert(dir == '+' || dir == '-');
    switch(dir) {
    case '+' : return_current = (is_ori ? true : false); break;
    case '-' : return_current = (is_ori ? false : true); break;
    default: assert(false);
    }
  }

  int getLen() const { return len; }

  int opAt(int pos) const {
    assert(pos >= 0 && pos < len);
    return bam_cigar_op(return_current ? cigar[pos] : cigar[len - pos - 1]);
  }

  char opchrAt(int pos) const { 
    assert(pos >= 0 && pos < len);
    return bam_cigar_opchr(return_current ? cigar[pos] : cigar[len - pos - 1]);
  }

  int oplenAt(int pos) const {
    assert(pos >= 0 && pos < len);
    return bam_cigar_oplen(return_current ? cigar[pos] : cigar[len - pos - 1]);
  }

  // toString will reset dir
  std::string toString(char dir = '+') {
    setDir(dir);
    std::ostringstream strout;
    for (int i = 0; i < len; ++i) strout<< opAt(i)<< oplenAt(i);
    return strout.str();
  }

private:
  const uint32_t *cigar;
  bool is_ori; // if the stored cigar is in original read sequence's order
  bool return_current; // if we can return the current cigar or return the reverse
  int len; // len, cigar string length
};

#endif
