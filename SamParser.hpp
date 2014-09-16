#ifndef SAMPARSER_H_
#define SAMPARSER_H_

#include "sam/sam.h"

#include "BamAlignment.hpp"
#include "AlignmentGroup.hpp"

class SamParser {
public:
  SamParser(char inpType, const char* inpF, const char* aux = NULL);
  ~SamParser();

  const bam_header_t* getHeader() const { 
    return header;
  }

  bool next(BamAlignment& b) { return b.read(sam_in); }

  bool next(AlignmentGroup& ag) { return ag.read(sam_in); }
 
private:
  samfile_t *sam_in;
  bam_header_t *header;
};

#endif

