#include<cctype>
#include<cstdio>
#include<cstdlib>
#include<cstring>
#include<map>
#include<vector>
#include<string>
#include<algorithm>
#include<iostream>

#include<stdint.h>
#include "sam/sam.h"
#include "sam/bam.h"

#include "utils.h"
#include "Transcript.hpp"
#include "Transcripts.hpp"
#include "TransBamAlignment.hpp"

using namespace std;


char tiF[STRLEN] = "human_hg19_0.ti";

Transcripts transcripts;
bam1_t *b;

void setName(bam1_t* b, const char* name) {
  int len = strlen(name);
  b->core.l_qname = len + 1;
  b->data_len = b->core.l_qname;
  b->data = (uint8_t*)calloc(b->data_len, sizeof(char));
  memcpy(b->data, name, b->core.l_qname);
}

void setCIGAR(bam1_t* b, const char* cigar) {
  vector<uint32_t> cigarStr;
  b->core.n_cigar = 0;
  int op, oplen;

  int i = 0, len = strlen(cigar);
  while (i < len) {
    switch(cigar[i]) {
    case 'M' : op = BAM_CMATCH; break;
    case 'I' : op = BAM_CINS; break;
    case 'D' : op = BAM_CDEL; break;
    default : assert(false);
    }
    ++i; oplen = 0;
    while (i < len && isdigit(cigar[i])) oplen = oplen * 10 + (cigar[i++] - '0');
    cigarStr.push_back(oplen << BAM_CIGAR_SHIFT | op);
    b->core.n_cigar++;
  }

  b->data_len += b->core.n_cigar * 4;
  b->data = (uint8_t*)realloc(b->data, b->data_len);
  uint32_t *p = bam1_cigar(b);
  for (int i = 0; i < b->core.n_cigar; i++) 
    p[i] = cigarStr[i];
}

void setSEQandQUAL(bam1_t* b, const char* seq, const char* qual) {
  b->core.l_qseq = strlen(seq);
  assert(b->core.l_qseq == strlen(qual));

  b->data_len += (b->core.l_qseq + 1) / 2 + b->core.l_qseq;
  b->data = (uint8_t*)realloc(b->data, b->data_len);
  uint8_t *s = bam1_seq(b);
  memset(s, 0, sizeof(uint8_t) * ((b->core.l_qseq + 1) / 2 + b->core.l_qseq));

  // Set seq
  for (int i = 0; i < b->core.l_qseq; i++) {
    uint8_t base = 0;
    switch(seq[i]) {
    case 'A' : base = 1; break;
    case 'C' : base = 2; break;
    case 'G' : base = 4; break;
    case 'T' : base = 8; break;
    case 'N' : base = 15; break;
    default: assert(false);
    }
    bam1_seq_seti(s, i, base);
  }

  // Set qual
  uint8_t *q = bam1_qual(b);
  for (int i = 0; i < b->core.l_qseq; i++) 
    q[i] = qual[i] - 33;
}

void printCIGAR(bam1_t* b) {
  uint32_t *c = bam1_cigar(b);
  for (int i = 0; i < b->core.n_cigar; i++) 
    printf("%c%d", bam_cigar_opchr(c[i]), bam_cigar_oplen(c[i]));
  printf("\n");
}

void printSEQ(bam1_t* b) {
  uint8_t *s = bam1_seq(b);
  for (int i = 0; i < b->core.l_qseq; i++) {
    char base;
    switch(bam1_seqi(s, i)) {
    case 1 : base = 'A'; break;
    case 2 : base = 'C'; break;
    case 4 : base = 'G'; break;
    case 8 : base = 'T'; break;
    case 15 : base = 'N'; break;
    default : assert(false);
    }
    printf("%c", base);
  }
  printf("\n");
}

void printQUAL(bam1_t* b) {
  uint8_t *q = bam1_qual(b);
  for (int i = 0; i < b->core.l_qseq; i++) 
    printf("%c", q[i] + 33);
  printf("\n");
}

int main() {

  transcripts.readFrom(tiF);

  b = bam_init1();
  b->core.pos = 1074;
  b->core.tid = 1;

  b->core.flag |= 0x0001;
  b->core.flag |= 0x0002;
  b->core.flag |= 0x0010;
  b->core.flag |= 0x0040;

  b->core.qual = 255;

  setName(b, "read1");
  setCIGAR(b, "I5M5I1D3M1");
  setSEQandQUAL(b, "AACGTATGCCTG", "#AAAAA!BBBBB");
  char MDstr[] = "2A2^TAT1";
  bam_aux_append(b, "MD", 'Z', strlen(MDstr) + 1, (uint8_t*)MDstr);
  printf("data_len = %d\n", b->data_len);
  float frac = 0.95;
  bam_aux_append(b, "ZW", 'f', bam_aux_type2size('f'), (uint8_t*)&frac);

  TransBamAlignment tba;

  // test if set bam works
  printf("Query Name = %s\n", bam1_qname(b));
  printf("CIGAR = "); printCIGAR(b);
  printf("SEQ = "); printSEQ(b);
  printf("QUAL = "); printQUAL(b);

  uint8_t *tag = NULL;
  tag = bam_aux_get(b, "MD");
  printf("MDstr = %s\n", bam_aux2Z(tag));
  tag = bam_aux_get(b, "ZW");
  printf("Frac = %f\n", bam_aux2f(tag));
  printf("\n\n");
  // Test is finished!

  // Test case 1
  tba.tr2chr(b, transcripts.getTranscriptAt(2));
  printf("Test 1:\n");
  printf("POS = %d\n", b->core.pos);
  printf("CIGAR = "); printCIGAR(b);
  printf("\n");
  bam_destroy1(b);

  // Test case 2;

  b = bam_init1();
  b->core.tid = 1;
  b->core.pos = 10;
  setName(b, "read2");
  setCIGAR(b, "M2130");
  
  char r2seq[2131], r2qual[2131];
  for (int i = 0; i < 2130; i++) {
    r2seq[i] = 'A';
    r2qual[i] = '#';
  }
  r2seq[2130] = r2qual[2130] = 0;
  
  setSEQandQUAL(b, r2seq, r2qual);

  tba.tr2chr(b, transcripts.getTranscriptAt(2));
  printf("Test 2:\n");
  printf("POS = %d\n", b->core.pos);
  printf("CIGAR = "); printCIGAR(b);
  printf("\n");
  bam_destroy1(b);

  // Test case 3;

  b = bam_init1();
  b->core.tid = 0;
  b->core.pos = 1535;
  
  setName(b, "read3");
  setCIGAR(b, "M5D3I2M230D5I4");
  char r3seq[242], r3qual[242];
  for (int i = 0; i < 241; i++) {
    r3seq[i] = 'C'; r3qual[i] = 'B';
  }
  r3seq[241] = r3qual[241] = 0;
  setSEQandQUAL(b, r3seq, r3qual);
  
  char MDstr3[] = "5^ACC110G119^TCGT";
  bam_aux_append(b, "MD", 'Z', strlen(MDstr3) + 1, (uint8_t*)MDstr3);

  tba.tr2chr(b, transcripts.getTranscriptAt(1));
  tba.modifyTags('-', b);
  printf("Test 3:\n");
  printf("POS = %d\n", b->core.pos);
  printf("CIGAR = "); printCIGAR(b);
  tag = bam_aux_get(b, "MD");
  printf("MDstr = %s\n", bam_aux2Z(tag));
  tag = bam_aux_get(b, "XS");
  printf("XS:A = %c\n", bam_aux2A(tag));
  printf("\n");
  bam_destroy1(b);

  // Test 4
  b = bam_init1();
  b->core.tid = 0;
  b->core.pos = 965; 
  
  setName(b, "read4");
  setCIGAR(b, "M10");
  setSEQandQUAL(b, "ACCGTATACC", "####ZZZZZZ");

  tba.tr2chr(b, transcripts.getTranscriptAt(1));
  tba.flipSeq(bam1_seq(b), b->core.l_qseq);
  tba.flipQual(bam1_qual(b), b->core.l_qseq);
  printf("Test 4:\n");
  printf("POS = %d\n", b->core.pos);
  printf("CIGAR = "); printCIGAR(b);
  printf("SEQ = "); printSEQ(b);
  printf("QUAL = "); printQUAL(b);
  printf("\n");
  bam_destroy1(b);

  return 0;
}
