#include<cstdio>
#include<cstring>
#include<cstdlib>
#include<cassert>

#include "sam/sam.h"
#include "sam/sam_header.h"

#include "my_assert.h"
#include "SamParser.hpp"

SamParser::SamParser(char inpType, const char* inpF, const char* aux) {
  switch(inpType) {
  case 'b': sam_in = samopen(inpF, "rb", aux); break;
  case 's': sam_in = samopen(inpF, "r", aux); break;
  default: assert(false);
  }

  general_assert(sam_in != 0, "Cannot open " + cstrtos(inpF) + "! It may not exist.");
  header = sam_in->header;
  general_assert(header != 0, "Fail to parse the header!");
}

SamParser::~SamParser() {
  samclose(sam_in);
}

const char* SamParser::getProgramID() {
  const char *key;
  if (header->dict == 0) header->dict = sam_header_parse2(header->text);                                                                                                                                 

  general_assert(sam_header2key(header->dict, "PG","ID", &key) != NULL, "Input SAM/BAM file does not contain a program ID!");
  return key;
}
