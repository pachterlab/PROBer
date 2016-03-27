#include<cstdio>
#include<cstring>
#include<cstdlib>
#include<string>
#include<sstream>

#include<stdint.h>
#include "htslib/sam.h"

#include "my_assert.h"
#include "BamWriter.hpp"

// If header == NULL, just create an empty header with one target (to avoid free non null pointer due to calloc) and paste the program_id line
// if program_id == NULL, use the header passed
BamWriter::BamWriter(const char* outF, const bam_hdr_t* header, const char* program_id) {
  if (header != NULL) {
    this->header = program_id == NULL ? bam_hdr_dup(header) : header_duplicate_without_text(header);
  }
  else {
    this->header = bam_hdr_init();
    this->header->n_targets = 1;
    this->header->target_len = new uint32_t[1];
    this->header->target_name = new char*[1];
    this->header->target_len[0] = 0;
    this->header->target_name[0] = new char[1];
    this->header->target_name[0][0] = 0;
  }

  if (program_id != NULL) {
    std::ostringstream strout;
    strout<< "@HD\tVN:1.4\tSO:unknown\n@PG\tID:"<< program_id<< std::endl;
    header_append_new_text(this->header, strout.str());
  }
  
  bam_out = sam_open(outF, "wb");
  general_assert(bam_out != 0, "Cannot write to " + cstrtos(outF) + "!");
  sam_hdr_write(bam_out, this->header);
}

BamWriter::~BamWriter() {
  bam_hdr_destroy(header);
  sam_close(bam_out);
}

bam_hdr_t* BamWriter::header_duplicate_without_text(const bam_hdr_t *ori_h) {
  bam_hdr_t *h = bam_hdr_init();
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

void BamWriter::header_append_new_text(bam_hdr_t* header, const std::string& new_text) {
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
