#ifndef HIT_H_
#define HIT_H_

struct Hit {
  int sid; // where the read aligns to, 0 means noise transcript 
  double prob; // the probability or fraction
};

#endif
