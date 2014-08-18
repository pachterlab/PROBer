#ifndef SEQSTRING_H_
#define SEQSTRING_H_

#include<cassert>
#include<string>
#include<sstream>

#include<stdint.h>
#include "sam/bam.h"

// An iterator for BAM sequence strings
class SEQstring {
public:
  SEQstring() : seq(NULL), dir(0), len(0) {}
  
  void setUp(uint8_t *seq, int len, char dir = '+') { 
    this->seq = seq;
    this->len = len;

    assert(dir == '+' || dir == '-');
    this->dir = dir;
  }

  void setDir(char dir) { 
    assert(dir == '+' || dir == '-');
    this->dir = dir;
  }

  int getLen() const { return len; }

  char baseAt(int pos) const {
    assert(pos >= 0 && pos < len);
    return (dir == '+' ? decode[bam1_seqi(seq, pos)] : decode_r[bam1_seqi(seq, len - pos - 1)]);
  }

  int baseCodeAt(int pos) const {
    assert(pos >= 0 && pos < len);
    return (dir == '+' ? codes[bam1_seqi(seq, pos)] : rcodes[bam1_seqi(seq, len - pos - 1)]);
  }

  // toString will reset dir
  std::string toString(char = '+');

private:
  uint8_t *seq;
  char dir;
  int len; // len, sequence length

  static const char decode[17], decode_r[17];
  static const int codes[16], rcodes[16];
};

#endif
