#ifndef BAMWRITER_H_
#define BAMWRITER_H_

#include<string>

#include "sam/bam.h"
#include "sam/sam.h"

#include "BamAlignment.hpp"
#include "AlignmentGroup.hpp"

class BamWriter {
public:
  BamWriter(const char* outF, const bam_header_t* header, const char* program_id);
  ~BamWriter();

  // choice: 0, do nothing; 1, delete read sequence and qual score; 2, add read sequence and qual score  
  bool write(BamAlignment& b, int choice = 0) { return b.write(bam_out, choice); }
  
  // choice: 0, do nothing; 1, delete read sequence and qual score; 2, add read sequence and qual score
  bool write(AlignmentGroup& ag, int choice = 0) { return ag.write(bam_out, choice); } 

private:
  samfile_t *bam_out;

  bam_header_t* header_duplicate_without_text(const bam_header_t* ori_h);
  void header_append_new_text(bam_header_t* header, const std::string& new_text);
};

#endif
