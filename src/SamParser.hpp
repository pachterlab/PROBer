#ifndef SAMPARSER_H_
#define SAMPARSER_H_

#include<string>

#include "htslib/sam.h"

#include "BamAlignment.hpp"
#include "AlignmentGroup.hpp"

class SamParser {
public:
  SamParser(const char* inpF, bam_hdr_t* input_header = NULL);
  ~SamParser();

  const bam_hdr_t* getHeader() const { return header; }

  bam_hdr_t* pass_header() { delete_header = false; return header; } // pass the header to an outside variable, then it is the outside variable's responsibility to delete header

  const char* getProgramID(); // scan header to look up program ID, slow

  bool next(BamAlignment& b) { return b.read(sam_in, header); }

  bool next(AlignmentGroup& ag) { return ag.read(sam_in, header); }
  
private:
  samFile* sam_in;
  bam_hdr_t* header;

  char program_id[1005];
  bool delete_header;
};

#endif

