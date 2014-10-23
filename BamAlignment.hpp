#ifndef BAMALIGNMENT_H_
#define BAMALIGNMENT_H_

#include<cmath>
#include<cassert>
#include<string>
#include<algorithm>

#include<stdint.h>
#include "sam/bam.h"
#include "sam/sam.h"

#include "my_assert.h"

#include "CIGARstring.hpp"
#include "SEQstring.hpp"
#include "QUALstring.hpp"

class BamAlignment {
public:
  BamAlignment() : is_paired(false), is_aligned(-1), b(NULL), b2(NULL) {
  }

  BamAlignment(const BamAlignment& o) : b(NULL), b2(NULL) {
    is_paired = o.is_paired;
    is_aligned = o.is_aligned;
    if (o.b != NULL) b = bam_dup1(o.b);
    if (o.b2 != NULL) b2 = bam_dup1(o.b2);
  }

  ~BamAlignment() {
    if (b != NULL) bam_destroy1(b);
    if (b2 != NULL) bam_destroy1(b2);
  }
  
  bool read(samfile_t*, BamAlignment* = NULL);
  bool write(samfile_t*, int = 0, BamAlignment* = NULL);
  
  // overall stats
  
  bool isPaired() const { return is_paired; }
  
  // -1: nothing is loaded, 0: not aligned, 1: fully aligned, 2: first mate is aligned (including SE reads), 3: second mate is aligned. 
  int isAligned() const { return is_aligned; }
  
  // mate = 0, the name of the read; 1, mate 1; 2, mate 2
  const char* getName(int mate = 0) const { 
    assert((mate == 0) || (mate == 1) || (mate == 2 && is_paired));
    return (char*)bam1_qname(mate < 2 ? b : b2); 
  }
  
  int getInsertSize() const {
    assert(is_aligned == 1);
    return abs(b->core.isize); 
  }
  
  int getTid() const { 
    assert(is_aligned > 0); 
    return (is_aligned != 3 ? b->core.tid : b2->core.tid);
  }
  
  char getDir() const { 
    assert(is_aligned > 0);
    if (is_aligned != 3) return !(b->core.flag & 0x0010) ? '+' : '-';
    return (b2->core.flag & 0x0010) ? '+' : '-';
  }

  /*
    @param     fragment_length     The average fragment length, 0 means no fragment length is provided
    @return    If fragment_length == 0, return the leftmost position of two mates. Otherwise, return the leftmost position calculated with fragment length.
               If the calculated position < 0, set it to 0
   */
  int getLeftMostPos(int fragment_length = 0) const {
    assert(is_aligned > 0);
    if (is_aligned == 1) return std::min(b->core.pos, b2->core.pos);
    if (fragment_length <= 0) return (is_aligned == 2 ? b->core.pos : b2->core.pos);

    if (is_aligned == 2) {
      return ((b->core.flag & 0x0010) == 0) ? b->core.pos : std::max(0, (int)bam_calend(&(b->core), bam1_cigar(b)) - fragment_length);
    }
    else {
      return ((b2->core.flag & 0x0010) == 0) ? b2->core.pos : std::max(0, (int)bam_calend(&(b2->core), bam1_cigar(b2)) - fragment_length);
    }
  }
  
  // if length = 2k, midpoint is k - 1
  int getMidPos(int fragment_length = 0) const {
    assert(is_aligned > 0);
    return (is_aligned == 1 ? getLeftMostPos() + (getInsertSize() - 1) / 2 : getLeftMostPos(fragment_length) + (fragment_length - 1) / 2);
  }
  
  int getMapQ() const { return (is_aligned != 3 ? b->core.qual : b2->core.qual); }
  
  void setMapQ(int MapQ) { 
    b->core.qual = ((b->core.flag & 0x0004) == 0 ? MapQ : 255);
    if (is_paired) b2->core.qual = ((b2->core.flag & 0x004) == 0 ? MapQ : 255);
  }
  
  // for mates
  
  // position from the strand where the mate comes from
  int getDirPos(int mate, const bam_header_t *header) const { 
    assert(is_aligned > 0);
    switch(mate) {
    case 1: return (b->core.flag & 0x0010) == 0 ? b->core.pos : header->target_len[b->core.tid] - b->core.pos - 1;
    case 2: 
      assert(is_paired);
      return (b2->core.flag & 0x0010) == 0 ? b2->core.pos : header->target_len[b2->core.tid] - b2->core.pos - 1; 
    default: assert(false);
    }
  }
  
  // length of the query sequence
  int getSeqLength(int mate = 1) const { 
    assert(mate == 1 || (is_paired && mate == 2));
    assert(((mate == 1) && (b->core.l_qseq > 0)) || ((mate == 2) && (b2->core.l_qseq > 0)));
    return mate == 1 ? b->core.l_qseq : b2->core.l_qseq;
  } 
  
  bool getCIGAR(CIGARstring& ci, int mate = 1) {
    assert(mate == 1 || (is_paired && mate == 2));
    if (mate == 1) ci.setUp(bam1_cigar(b), b->core.n_cigar, getSign(b));
    else ci.setUp(bam1_cigar(b2), b2->core.n_cigar, getSign(b2));
    return true;
  }

  bool getSEQ(SEQstring& si, int mate = 1) {
    assert(mate == 1 || (is_paired && mate == 2));
    assert(((mate == 1) && (b->core.l_qseq > 0)) || ((mate == 2) && (b2->core.l_qseq > 0)));
    if (mate == 1) si.setUp(bam1_seq(b), b->core.l_qseq, getSign(b));
    else si.setUp(bam1_seq(b2), b2->core.l_qseq, getSign(b2));
    return true;
  }
  
  bool getQUAL(QUALstring& qi, int mate = 1) {
    assert(mate == 1 || (is_paired && mate == 2));
    assert(((mate == 1) && (b->core.l_qseq > 0)) || ((mate == 2) && (b2->core.l_qseq > 0)));
    if (mate == 1) qi.setUp(bam1_qual(b), b->core.l_qseq, getSign(b));
    else qi.setUp(bam1_qual(b2), b2->core.l_qseq, getSign(b2));
    return true;
  }

  // optional fields

  void removeTag(const char tag[2]) {
    uint8_t *p_tag = NULL;

    while ((p_tag = bam_aux_get(b, tag)) != NULL) bam_aux_del(b, p_tag);
    if (is_paired) 
      while ((p_tag = bam_aux_get(b2, tag)) != NULL) bam_aux_del(b2, p_tag);
  }

  /*
    @func   this function append a ZF:A:! field to indicate this alignment is filtered out if it is not marked before
   */
  void markAsFiltered() {
    char c = '!';
    uint8_t *p_tag = NULL;

    p_tag = bam_aux_get(b, "ZF");
    if (p_tag == NULL) bam_aux_append(b, "ZF", 'A', bam_aux_type2size('A'), (uint8_t*)&c);
    if (is_paired) {
      p_tag = bam_aux_get(b2, "ZF");
      if (p_tag == NULL) bam_aux_append(b2, "ZF", 'A', bam_aux_type2size('A'), (uint8_t*)&c);
    }
  }

  bool findTag(const char tag[2], uint8_t*& p, char& type, int mate = 1) {
    assert(mate == 1 || (is_paired && mate == 2));
    p = bam_aux_get((mate == 1 ? b : b2), tag);
    if (p == NULL) return false;
    type = *p;
    if (type == 'c' || type == 'C' || type == 's' || type == 'S' || type == 'I') type = 'i';
    assert(type == 'A' || type == 'i' || type == 'f' || type == 'd' || type == 'Z' || type == 'H' || type == 'B');
    return true;
  }

  // no check, must guarantee the type is consistent
  char tag2A(uint8_t* p) { return bam_aux2A(p); }
  int tag2i(uint8_t* p) { return bam_aux2i(p); }
  float tag2f(uint8_t* p) { return bam_aux2f(p); }
  double tag2d(uint8_t* p) { return bam_aux2d(p); }
  char* tag2Z(uint8_t* p) { return bam_aux2Z(p); }
  char* tag2H(uint8_t* p) { return bam_aux2Z(p); }
  char* tag2B(uint8_t* p) { return (char*)(p + 1); }

  // If no ZW field, return -1.0
  float getFrac() {
    uint8_t *p_tag = bam_aux_get((is_aligned != 3 ? b : b2), "ZW");
    return p_tag != NULL ? bam_aux2f(p_tag) : -1.0;
  }
  
  void setFrac(float frac) {
    uint8_t *p_tag = NULL;

    if (!(b->core.flag & 0x0004)) {
      p_tag = bam_aux_get(b, "ZW");
      if (p_tag != NULL) {
	memcpy(p_tag + 1, (uint8_t*)&(frac), bam_aux_type2size('f'));
      } else {
	bam_aux_append(b, "ZW", 'f', bam_aux_type2size('f'), (uint8_t*)&frac);
      }
    }
    
    if (is_paired && !(b2->core.flag & 0x0004)) {
      p_tag = bam_aux_get(b2, "ZW");
      if (p_tag != NULL) {
	memcpy(p_tag + 1, (uint8_t*)&(frac), bam_aux_type2size('f'));
      } else {
	bam_aux_append(b2, "ZW", 'f', bam_aux_type2size('f'), (uint8_t*)&frac);
      }
    }

    setMapQ(frac2MapQ(frac));
  }

  void writeToBED(FILE*, const bam_header_t*);
  void writeToTagAlign(FILE*, const bam_header_t*);

protected:

  bool is_paired;
  char is_aligned; // 0, unalignable; 1, both mates are aligned; 2, only first mate is aligned; 3, only second mate is aligned
  
  bam1_t *b, *b2;
  
  uint8_t frac2MapQ(float val) {
    float err = 1.0 - val;
    if (err <= 1e-10) return 100;
    return (uint8_t)(-10 * log10(err) + .5); // round it
  }

  // Caution: this function may change b->data's adddress!  
  void expand_data_size(bam1_t *b) {
    if (b->m_data < b->data_len) {
      b->m_data = b->data_len;
      kroundup32(b->m_data);
      b->data = (uint8_t*)realloc(b->data, b->m_data);
    }
  }  
  
  char getSign(bam1_t *b) {
    return (((b->core.flag & 0x0004) || !(b->core.flag & 0x0010)) ? '+' : '-');
  }
};

#endif
