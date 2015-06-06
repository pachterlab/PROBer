#ifndef BAMCONVERTER_H_
#define BAMCONVERTER_H_

#include<string>
#include<map>

#include "sam/bam.h"
#include "sam/sam.h"

#include "Transcripts.hpp"
#include "TransAlignmentGroup.hpp"

class BamConverter {
public:
  BamConverter(const char*, const char*, const char*, Transcripts&);
  ~BamConverter();

  void process();
private:
  samfile_t *in, *out;
  Transcripts& transcripts;

  std::map<std::string, int> refmap;
  TransAlignmentGroup *tag;

  void header_append_new_text(bam_header_t*, const std::string&);
};

#endif
