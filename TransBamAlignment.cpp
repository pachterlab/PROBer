#include<cstdio>
#include<cstring>
#include<cstdlib>
#include<vector>

#include<stdint.h>
#include "sam/bam.h"

#include "utils.h"
#include "my_assert.h"
#include "Transcript.hpp"
#include "CIGARstring.hpp"
#include "BamAlignment.hpp"
#include "TransBamAlignment.hpp"

void TransBamAlignment::convert(int cid, const Transcript& transcript) {
  if (is_aligned == 0) return;

  if (!(b->core.flag & 0x0004)) convertBam1_t(b, cid, transcript);

  if (is_paired && !(b2->core.flag & 0x0004)) convertBam1_t(b2, cid, transcript);

  if (is_aligned == 1) {
    b->core.mtid = b2->core.tid;
    b->core.mpos = b2->core.pos;
    b2->core.mtid = b->core.tid;
    b2->core.mpos = b->core.pos;
  }
}

void TransBamAlignment::convertBam1_t(bam1_t *b, int cid, const Transcript& transcript) {
  char strand = transcript.getStrand();

  b->core.tid = cid;
  if (strand == '-') {
    b->core.flag ^= 0x0010;
    if (is_aligned == 1) {
      b->core.flag ^=0x0020;
      b->core.isize = -b->core.isize;
    }
  }

  tr2chr(b, transcript);

  if (strand == '-') {
    flipSeq(bam1_seq(b), b->core.l_qseq);
    flipQual(bam1_qual(b), b->core.l_qseq);
  }
  
  b->core.bin = bam_reg2bin(b->core.pos, bam_calend(&(b->core), bam1_cigar(b)));
  modifyTags(strand, b);
}

//convert transcript coordinate to chromosome coordinate and generate CIGAR string
void TransBamAlignment::tr2chr(bam1_t* b, const Transcript& transcript) {
  char strand = transcript.getStrand();
  int length = transcript.getLength();
  const std::vector<Interval>& structure = transcript.getStructure();
  int s = structure.size();
  assert(s > 0);

  std::vector<uint32_t> cigarStr;
  cigarStr.clear();
  
  if (b->core.pos >= length) { // Read aligned to the polyA tail totally
    if (strand == '+') b->core.pos = structure[s - 1].end; // set the start position be the last chromosome position + 1, 0-based
    else b->core.pos = structure[0].start - 2; // set the start position be the first chromosome position - 1, 0-based
    cigarStr.push_back(b->core.l_qseq << BAM_CIGAR_SHIFT | BAM_CINS);
  }
  else {
    CIGARstring cigar;
    cigar.setUp(bam1_cigar(b), b->core.n_cigar);
    cigar.setDir(strand);
    int sp = (strand == '+' ? b->core.pos : length - bam_calend(&(b->core), bam1_cigar(b))); // 0-based, bam_calend gives the right most position + 1
    int len = cigar.getLen();
    assert(len > 0);

    int curlen = 0;
    int op, oplen = 0, i = -1; // oplen, the left length of current operation; i, the position of the current cigar operation  

    int minlen;

    // load the first cigar operation
    ++i; op = cigar.opAt(i); oplen = cigar.oplenAt(i);
    general_assert(op == BAM_CINS || op == BAM_CDEL || op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF, "Detected cigar operation other than I/D/M/=/X in a transcript alignment!"); 

    // If the beginning of the read aligns to the poly(A) tail
    while (sp < 0 && i < len) {
      if (op == BAM_CINS) { 
	curlen += oplen; oplen = 0;
      }
      else {
	minlen = std::min(-sp, oplen);
	sp += minlen; oplen -= minlen;
	if (op != BAM_CDEL) curlen += minlen;
      }

      if (oplen == 0) {
	++i; if (i >= len) continue;
	op = cigar.opAt(i); oplen = cigar.oplenAt(i); 
	general_assert(op == BAM_CINS || op == BAM_CDEL || op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF, "Detected cigar operation other than I/D/M/=/X in a transcript alignment!"); 
      }
    }

    // Remove leading insertions
    while (i < len && op == BAM_CINS) {
      curlen += oplen; 
      ++i; if (i >= len) continue;
      op = cigar.opAt(i); oplen = cigar.oplenAt(i); 
      general_assert(op == BAM_CINS || op == BAM_CDEL || op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF, "Detected cigar operation other than I/D/M/=/X in a transcript alignment!");
    }

    assert(i < len);

    if (curlen > 0) cigarStr.push_back(curlen << BAM_CIGAR_SHIFT | BAM_CINS);

    // Determine the leftmost genomic coordinate of the read
    int strup = 0, strulen = structure[0].end - structure[0].start + 1;

    while (strulen <= sp && strup + 1 < s) { ++strup; strulen += structure[strup].end - structure[strup].start + 1; }
    assert(strulen > sp); 

    b->core.pos = structure[strup].end - (strulen - sp);
    
    // Determine cigar strings
    while (i < len && strup < s) {
      if (op == BAM_CINS) {
	cigarStr.push_back(oplen << BAM_CIGAR_SHIFT | op);
	oplen = 0;
      }
      else {
	assert(sp <= strulen); // guard test, comment out before release
	if (sp == strulen) {
	  ++strup; if (strup >= s) continue;
	  curlen = structure[strup].start - structure[strup - 1].end - 1;
	  cigarStr.push_back(curlen << BAM_CIGAR_SHIFT | BAM_CREF_SKIP);
	  strulen += structure[strup].end - structure[strup].start + 1;
	}
	minlen = std::min(strulen - sp, oplen);
	cigarStr.push_back(minlen << BAM_CIGAR_SHIFT | op);
	sp += minlen;
	oplen -= minlen;
      }
      
      if (oplen == 0) {
	++i; if (i >= len) continue;
	op = cigar.opAt(i); oplen = cigar.oplenAt(i); 
	general_assert(op == BAM_CINS || op == BAM_CDEL || op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF, "Detected cigar operation other than I/D/M/=/X in a transcript alignment!");
      }
    }

    // the end of the read aligns to the poly(A) tail
    curlen = 0;
    while (i < len) {
      if (op != BAM_CDEL) curlen += oplen;
      ++i; if (i >= len) continue;
      op = cigar.opAt(i); oplen = cigar.oplenAt(i); 
      general_assert(op == BAM_CINS || op == BAM_CDEL || op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF, "Detected cigar operation other than I/D/M/=/X in a transcript alignment!");
    }

    assert(cigarStr.size() > 0);

    if (curlen > 0) {
      if (bam_cigar_op(cigarStr.back()) == BAM_CINS) {
	curlen += bam_cigar_oplen(cigarStr.back());
	cigarStr.pop_back();
      }
      cigarStr.push_back(curlen << BAM_CIGAR_SHIFT | BAM_CINS);
    }
  }

  // write down the new CIGAR string
  b->data_len += ((int)cigarStr.size() - b->core.n_cigar) * 4;
  expand_data_size(b); // Caution: this function may change b->data's adddress!
  uint8_t *seq = bam1_seq(b);
  int num = (b->core.l_qseq + 1) / 2 + b->core.l_qseq + b->l_aux;
  b->core.n_cigar = cigarStr.size();
  memmove(bam1_seq(b), seq, num);
  uint32_t *c = bam1_cigar(b);
  for (int i = 0; i < b->core.n_cigar; i++) c[i] = cigarStr[i];
}

void TransBamAlignment::flipSeq(uint8_t* s, int32_t l_qseq) {
  uint8_t code, base;
  uint8_t *seq = new uint8_t[(l_qseq + 1) / 2];
  int32_t len = 0;

  code = 0; base = 0;
  for (int32_t i = 0; i < l_qseq; i++) {
    switch (bam1_seqi(s, l_qseq - i - 1)) {
    case 1: base = 8; break;
    case 2: base = 4; break;
    case 4: base = 2; break;
    case 8: base = 1; break;
    case 15: base = 15; break;
    default: assert(false);
    }
    code |=  base << ((1 - (i & 1)) << 2);
    if (i & 1) { seq[len++] = code; code = 0; }
  }
  if (l_qseq & 1) seq[len++] = code;
  memcpy(s, seq, len);
  delete[] seq;
}

void TransBamAlignment::flipQual(uint8_t* q, int32_t l_qseq) {
  int32_t mid = l_qseq / 2;
  uint8_t tmp;
  for (int32_t i = 0; i < mid; i++) {
    tmp = q[i]; q[i] = q[l_qseq - i - 1]; q[l_qseq - i - 1] = tmp;
  }
}

// MD field might not match the CIGAR field for the portion that aligns to the poly(A) tail. Any deletion at the poly(A) tail will be ommitted but still be reflected in the MD field
void TransBamAlignment::modifyTags(char strand, bam1_t* b) {
  uint8_t *s = NULL;
  
  if (strand == '-') {
    s = bam_aux_get(b, "MD");
    if ((s != NULL) && (*(s) == 'Z') && (bam_aux2Z(s) != NULL)) {
      char *mis = bam_aux2Z(s);
      int len = strlen(mis);
      char *tmp = new char[len];
      int cur_type = -1, fr = -1, type, base;
      for (int i = 0; i < len; i++) {
	type = (mis[i] >= '0' && mis[i] <= '9');
	if (cur_type != type) {
	  switch(cur_type) {
	  case 0:
	    base = len - 1;
	    if (mis[fr] == '^') { tmp[len - i] = mis[fr]; ++fr; ++base; }
	    for (int j = fr; j < i; j++) tmp[base - j] = ((mis[j] == 'A' || mis[j] == 'C' || mis[j] == 'G' || mis[j] == 'T') ? getOpp(mis[j]) : mis[j]);
	    break;
	  case 1: 
	    base = len - i - fr;
	    for (int j = fr; j < i; j++) tmp[base + j] = mis[j]; 
	    break; 
	  }
	  cur_type = type;
	  fr = i;
	}
      }
      switch(cur_type) {
      case 0:
	base = len - 1;
	if (mis[fr] == '^') { tmp[0] = mis[fr]; ++fr; ++base; }
	for (int j = fr; j < len; j++) tmp[base - j] = ((mis[j] == 'A' || mis[j] == 'C' || mis[j] == 'G' || mis[j] == 'T') ? getOpp(mis[j]) : mis[j]);
	break;
      case 1: 
	for (int j = fr; j < len; j++) tmp[j - fr] = mis[j]; 
	break; 
      }
      strncpy(mis, tmp, len);
      delete[] tmp;
    }
  }
  
  // append XS:A field if necessary
  s = bam_aux_get(b, "XS");
  if (s != NULL) bam_aux_del(b, s);
  bool hasN = false;
  uint32_t* p = bam1_cigar(b);
  for (int i = 0; i < (int)b->core.n_cigar; i++)
    if ((*(p + i) & BAM_CIGAR_MASK) == BAM_CREF_SKIP) { hasN = true; break; }
  if (hasN) bam_aux_append(b, "XS", 'A', 1, (uint8_t*)&strand);
}
