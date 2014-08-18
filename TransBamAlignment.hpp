#ifndef TRANSBAMALIGNMENT_H_
#define TRANSBAMALIGNMENT_H_

#include<stdint.h>
#include "sam/bam.h"

#include "Transcript.hpp"
#include "BamAlignment.hpp"

class TransBamAlignment : public BamAlignment {
public:
  TransBamAlignment() : BamAlignment() {}
  TransBamAlignment(const TransBamAlignment& o) : BamAlignment(o) {}

  int compare(const TransBamAlignment* o) const {
    int value;
    if (is_aligned != o->is_aligned) return getSign(is_aligned < o->is_aligned);
    if (is_aligned != 3) {
      value = compareBam1_t(b, o->b);
      if (value != 0) return value;
    }
    if (is_aligned != 2) {
      value = compareBam1_t(b2, o->b2);
      if (value != 0) return value;
    }
    return 0;
  }

  void convert(int, const Transcript&);

private:

  int getSign(bool value) const { return value ? -1 : 1; }

  int compareBam1_t(const bam1_t *b, const bam1_t *ob) const {
    int strand1, strand2;
    uint32_t *p1, *p2;
    
    if (b->core.tid != ob->core.tid) return getSign(b->core.tid < ob->core.tid);
    if (b->core.pos != ob->core.pos) return getSign(b->core.pos < ob->core.pos);
    strand1 = b->core.flag & 0x0010; strand2 = ob->core.flag & 0x0010;
    if (strand1 != strand2) return getSign(strand1 < strand2);
    if (b->core.n_cigar != ob->core.n_cigar) return getSign(b->core.n_cigar < ob->core.n_cigar);
    p1 = bam1_cigar(b); p2 = bam1_cigar(ob);
    for (int i = 0; i < (int)b->core.n_cigar; i++) {
      if (*p1 != *p2) return getSign(*p1 < *p2);
      ++p1; ++p2;
    }
    
    return 0;
  }

  void convertBam1_t(bam1_t*, int cid, const Transcript&);
  void tr2chr(bam1_t*, const Transcript&);
  void flipSeq(uint8_t*, int32_t);
  void flipQual(uint8_t*, int32_t);
  void modifyTags(char, bam1_t*);
  
};

#endif

