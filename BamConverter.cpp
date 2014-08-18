#include<cstdio>
#include<cstring>
#include<cstdlib>
#include<cassert>
#include<string>
#include<sstream>
#include<map>

#include "sam/bam.h"
#include "sam/sam.h"

#include "utils.h"
#include "my_assert.h"
#include "Transcripts.hpp"
#include "TransAlignmentGroup.hpp"
#include "BamConverter.hpp"

BamConverter::BamConverter(const char* inpF, const char* outF, const char* chr_list, Transcripts& transcripts) : transcripts(transcripts) {
  general_assert(transcripts.getType() == 0, "Genome information is not provided! RSEM cannot convert the transcript bam file!");
  
  in = samopen(inpF, "rb", NULL);
  assert(in != 0);
  
  transcripts.buildMappings(NULL, in->header->n_targets, in->header->target_name);
  
  bam_header_t *out_header = sam_header_read2(chr_list);
  
  refmap.clear();
  for (int i = 0; i < out_header->n_targets; i++) {
    refmap[out_header->target_name[i]] = i;
  }
  
  if (in->header->l_text > 0) {
    std::istringstream strin(in->header->text);
    std::ostringstream strout;
    std::string line;
    
    while (std::getline(strin, line)) {
      if (line.substr(0, 3) != "@SQ") strout<< line<< "\n";
    }

    header_append_new_text(out_header, strout.str());
  }

  header_append_new_text(out_header, "@CO\tThis BAM file is processed by rsem-tbam2gam to convert from transcript coordinates into genomic coordinates.\n");

  out = samopen(outF, "wb", out_header);
  assert(out != 0);
  
  bam_header_destroy(out_header);  

  tag = NULL;
}

BamConverter::~BamConverter() {
  samclose(in);
  samclose(out);
}

void BamConverter::process() {
  tag = new TransAlignmentGroup(transcripts, refmap);
  HIT_INT_TYPE cnt = 0;
  while (tag->read(in)) {
    tag->convert(out);
    ++cnt;
    if (cnt % 1000000 == 0) { printf("."); fflush(stdout); }
  }
  printf("\n");
  delete tag;
}

void BamConverter::header_append_new_text(bam_header_t* header, const std::string& new_text) {
  if (new_text == "") return;
  int len = new_text.length();
  int max_size = (header->text == NULL ? 0 : header->l_text + 1);
  kroundup32(max_size);
  if (max_size < int(header->l_text + len + 1)) {
    max_size = header->l_text + len + 1;
    kroundup32(max_size);
    header->text = (char*)realloc(header->text, max_size);
  }
  strcpy(header->text + header->l_text, new_text.c_str());
  header->l_text += len;
}
