#ifndef SAMPARSER_H_
#define SAMPARSER_H_

#include<string>

#include "htslib/sam.h"

#include "BamAlignment.hpp"
#include "AlignmentGroup.hpp"

class SamParser {
public:
  SamParser(const char* inpF);
  ~SamParser();

  const bam_hdr_t* getHeader() const { return header; }

  const char* getProgramID(); // scan header to look up program ID, slow

  bool next(BamAlignment& b) { return b.read(sam_in, header); }

  bool next(AlignmentGroup& ag) { return ag.read(sam_in, header); }
  
private:
  samFile* sam_in;
  bam_hdr_t* header;

  const char program_id[1005];
};

#endif

