#include<cstdio>
#include<cstring>
#include<cstdlib>
#include<cassert>

#include "htslib/sam.h"

#include "my_assert.h"
#include "SamParser.hpp"

SamParser::SamParser(const char* inpF) {
  sam_in = sam_open(inpF, "r");
  general_assert(sam_in != 0, "Cannot open " + cstrtos(inpF) + "! It may not exist.");

  header = sam_hdr_read(sam_in);
  general_assert(header != 0, "Fail to parse the header!");

  memset(program_id, 0, sizeof(program_id));
}

SamParser::~SamParser() {
  bam_hdr_destroy(header);
  sam_close(sam_in);
}

// This is an simple implementation, improve it later if necessary
const char* SamParser::getProgramID() {
  if (program_id[0]) return program_id;
  
  char *p = strstr(header->text, "@PG\t");
  assert(p != NULL);
  p += 4;

  char *fr = p;
  while (*p != '\n' && *p != '\0') {
    if (*p == '\t') {
      if (p - fr > 3 && !strncmp(fr, "ID:", 3)) {
	assert(p - fr - 3 <= 100);
	strncpy(program_id, fr + 3, p - fr - 3);
	return program_id;
      }
      fr = p + 1;
    }
    ++p;
  }
  if (p - fr > 3 && !strncmp(fr, "ID:", 3)) {
    assert(p - fr - 3 <= 100);
    strncpy(program_id, fr + 3, p - fr - 3);
    return program_id;
  }
  assert(false);
}
