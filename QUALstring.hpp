#ifndef QUALSTRING_H_
#define QUALSTRING_H_

#include<cassert>
#include<string>
#include<sstream>

#include<stdint.h>

// An iterator for BAM quality score strings
class QUALstring {
public:
  QUALstring() : qual(NULL), dir(0), len(0) {}
  
  void setUp(uint8_t *qual, int len, char dir = '+') { 
    this->qual = qual;
    this->len = len;

    assert(dir == '+' || dir == '-');
    this->dir = dir;
  }

  void setDir(char dir) { 
    assert(dir == '+' || dir == '-');
    this->dir = dir;
  }

  int getLen() const { return len; }

  // 33 is already deducted
  int qualAt(int pos) const {
    assert(pos >= 0 && pos < len);
    return qual[dir == '+' ? pos : len - pos - 1];
  }

  // toString will reset dir
  std::string toString(char dir = '+') {
    assert(dir == '+' || dir == '-');
    this->dir = dir;
    std::ostringstream strout;
    for (int i = 0; i < len; i++) strout<< char(qualAt(i) + 33);
    this->dir = 0;
    return strout.str();
  }

private:
  uint8_t *qual;
  char dir;
  int len; // len, sequence length
};

#endif
