#ifndef MDSTRING_H_
#define MDSTRING_H_

#include<cstdio>
#include<cctype>
#include<cassert>

// An iterator for MD string of a BAM alignment
class MDstring {
public:
  MDstring() : mdstr(NULL) {}

  void setUp(const char* mdstr) {
    this->mdstr = mdstr;
    counter = 0;
  }

  // return next reference base
  char next() {
    if (counter > 0) { --counter; return 0; } // MATCH
    if (isdigit(*mdstr)) {
      counter = *mdstr - '0', ++mdstr;
      while (isdigit(*mdstr)) counter = counter * 10 + (*mdstr - '0'), ++mdstr;
      if (counter > 0) { --counter; return 0; }
    }
    if (*mdstr == 0) return -1; // end of the MD string
    while (*mdstr == '^') ++mdstr;
    assert(isalpha(*mdstr));
    return *mdstr++;
  }
  
private:
  const char* mdstr;
  int counter;
};


#endif
