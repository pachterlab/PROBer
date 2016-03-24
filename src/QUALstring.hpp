#ifndef QUALSTRING_H_
#define QUALSTRING_H_

#include<cassert>
#include<string>
#include<sstream>

#include<stdint.h>

// An iterator for BAM quality score strings
class QUALstring {
public:
  QUALstring() : qual(NULL), is_ori(true), return_current(true), len(0) {}
  
  void setUp(uint8_t *qual, int len, bool is_ori) { 
    this->qual = qual;
    this->len = len;
    this->is_ori = is_ori;
    // Default, return the original quality score sequence
    return_current = (is_ori ? true : false);
  }

  char getDir() { return ((is_ori && return_current) || (!is_ori && !return_current)) ? '+' : '-'; }
  
  // '+' returns the original qual string; '-' returns the reverse string
  void setDir(char dir) { 
    assert(dir == '+' || dir == '-');
    switch(dir) {
    case '+' : return_current = (is_ori ? true : false); break;
    case '-' : return_current = (is_ori ? false : true); break;
    default: assert(false);
    }
  }
  
  int getLen() const { return len; }
  
  // 33 is already deducted
  int qualAt(int pos) const {
    assert(pos >= 0 && pos < len);
    return qual[return_current ? pos : len - pos - 1];
  }

  // default is the original quality score string
  std::string toString(char dir = '+') {
    setDir(dir);
    std::ostringstream strout;
    for (int i = 0; i < len; i++) strout<< char(qualAt(i) + 33);
    return strout.str();
  }

private:
  uint8_t *qual;
  bool is_ori; // if the stored qual is in read's original form
  bool return_current; // if we can return the current seq or return the reverse
  int len; // len, quality score sequence length
};

#endif
