#include<cstdio>
#include<cstring>
#include<cstdlib>
#include<string>
#include<sstream>

#include<stdint.h>
#include "sam/bam.h"
#include "sam/sam.h"

#include "my_assert.h"

#include "BamWriter.hpp"

// If header == NULL, just create an empty header and paste the program_id line
BamWriter::BamWriter(const char* outF, const bam_header_t* header, const char* program_id) {
  bam_header_t *out_header = (header != NULL ? header_duplicate_without_text(header) : bam_header_init());
  
  std::ostringstream strout;
  strout<< "@HD\tVN:1.4\tSO:unknown\n@PG\tID:"<< program_id<< std::endl;
  header_append_new_text(out_header, strout.str());

  bam_out = samopen(outF, "wb", out_header);
  general_assert(bam_out != 0, "Cannot write to " + cstrtos(outF) + "!");

  bam_header_destroy(out_header);
}

BamWriter::~BamWriter() {
  samclose(bam_out);
}

bam_header_t* BamWriter::header_duplicate_without_text(const bam_header_t *ori_h) {
  bam_header_t *h = bam_header_init();
  h->n_targets = ori_h->n_targets;
  h->target_len = new uint32_t[h->n_targets];
  h->target_name = new char*[h->n_targets];
  for (int i = 0; i < h->n_targets; ++i) {
    h->target_len[i] = ori_h->target_len[i];
    h->target_name[i] = new char[strlen(ori_h->target_name[i]) + 1];
    strcpy(h->target_name[i], ori_h->target_name[i]);
  }
  return h;
}

void BamWriter::header_append_new_text(bam_header_t* header, const std::string& new_text) {
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
