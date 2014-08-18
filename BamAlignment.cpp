#include<cmath>
#include<cstdio>
#include<cstring>
#include<cassert>
#include<string>
#include<algorithm>

#include<stdint.h>
#include "sam/bam.h"
#include "sam/sam.h"

#include "my_assert.h"
#include "SEQstring.hpp"
#include "BamAlignment.hpp"

bool BamAlignment::read(samfile_t *in, BamAlignment *o) {
  is_aligned = -1;

  if (b == NULL) b = bam_init1();
  if (samread(in, b) < 0) return false;

  is_paired = (b->core.flag & 0x0001) > 0;
  // read the second mate
  if (is_paired) { 
    if (b2 == NULL) b2 = bam_init1();

    general_assert(samread(in, b2) >= 0 && (b2->core.flag & 0x0001), "Fail to read the other mate for a paired-end alignment!");
    general_assert(((b->core.flag & 0x00C0) == 0x0040 && (b2->core.flag & 0x00C0) == 0x0080) || 
		   ((b->core.flag & 0x00C0) == 0x0080 && (b2->core.flag & 0x00C0) == 0x0040), 
		   "Cannot detect both mates of a paired-end alignment!");
    
    if (b->core.flag & 0x0080) { bam1_t* tmp = b; b = b2; b2 = tmp; }
  }

  // calculate is_aligned
  if (!is_paired) is_aligned = ((b->core.flag & 0x0004) ? 0 : 2); // SE, 0 or 2
  else if (b->core.flag & 0x0004) is_aligned = ((b2->core.flag & 0x0004) ? 0 : 3);
  else is_aligned = ((b2->core.flag & 0x0004) ? 2 : 1);
  
  if (b->core.l_qseq <= 0) b->core.l_qseq = o->getSeqLength(1);
  if (is_paired && b2->core.l_qseq <= 0) b2->core.l_qseq = o->getSeqLength(2);

  // May consider to comment the following two lines out for efficiency
  assert((is_aligned == 0 || is_aligned == 3) || b->core.l_qseq == bam_cigar2qlen(&b->core, bam1_cigar(b)));
  assert((is_aligned == 0 || is_aligned == 2) || b2->core.l_qseq == bam_cigar2qlen(&b2->core, bam1_cigar(b2)));

  return true;
}

// choice: 0, do nothing; 1, delete read sequence and qual score; 2, add read sequence and qual score. o, the alignment that contain read sequence and quality score information
bool BamAlignment::write(samfile_t *out, int choice, BamAlignment *o) {
  assert(is_aligned >= 0 && b != NULL && (!is_paired || b2 != NULL));

  switch(choice) {
  case 0: 
    if (b->core.l_qname == 1) b->core.l_qseq = 0;
    if (is_paired && (b2->core.l_qname == 1)) b2->core.l_qseq = 0;
    break;
  case 1: 
    if (b->core.l_qname > 1) {
      b->data[0] = 0;
      memmove(b->data + 1, bam1_cigar(b), b->core.n_cigar * 4);
      memmove(b->data + 1 + b->core.n_cigar * 4, bam1_aux(b), b->l_aux);
      b->data_len = 1 + b->core.n_cigar * 4 + b->l_aux;
      b->core.l_qname = 1;
      b->core.l_qseq = 0;
    }
    if (is_paired && (b2->core.l_qname > 1)) {
      b2->data[0] = 0;
      memmove(b2->data + 1, bam1_cigar(b2), b2->core.n_cigar * 4);
      memmove(b2->data + 1 + b2->core.n_cigar * 4, bam1_aux(b2), b2->l_aux);
      b2->data_len = 1 + b2->core.n_cigar * 4 + b2->l_aux;
      b2->core.l_qname = 1;
      b2->core.l_qseq = 0;
    }
    break;
  case 2:
    if (b->core.l_qname == 1) {
      b->core.l_qname = o->b->core.l_qname;
      b->core.l_qseq = o->b->core.l_qseq;
      b->data_len = b->core.l_qname + b->core.n_cigar * 4 + (b->core.l_qseq + 1) / 2 + b->core.l_qseq + b->l_aux;
      expand_data_size(b);
      memmove(bam1_aux(b), b->data + 1 + b->core.n_cigar * 4, b->l_aux); // move aux options
      memmove(bam1_cigar(b), b->data + 1, b->core.n_cigar * 4); // move cigar string
      memcpy(bam1_qname(b), bam1_qname(o->b), b->core.l_qname); // copy qname
      memcpy(bam1_seq(b), bam1_seq(o->b), (b->core.l_qseq + 1) / 2); // copy seq
      memcpy(bam1_qual(b), bam1_qual(o->b), b->core.l_qseq); // copy qual
    }
    if (is_paired && (b2->core.l_qname == 1)) {
      b2->core.l_qname = o->b2->core.l_qname;
      b2->core.l_qseq = o->b2->core.l_qseq;
      b2->data_len = b2->core.l_qname + b2->core.n_cigar * 4 + (b2->core.l_qseq + 1) / 2 + b2->core.l_qseq + b2->l_aux;
      expand_data_size(b2);
      memmove(bam1_aux(b2), b2->data + 1 + b2->core.n_cigar * 4, b2->l_aux); // move aux options
      memmove(bam1_cigar(b2), b2->data + 1, b2->core.n_cigar * 4); // move cigar string
      memcpy(bam1_qname(b2), bam1_qname(o->b2), b2->core.l_qname); // copy qname
      memcpy(bam1_seq(b2), bam1_seq(o->b2), (b2->core.l_qseq + 1) / 2); // copy seq
      memcpy(bam1_qual(b2), bam1_qual(o->b2), b2->core.l_qseq); // copy qual     
    }
    break;
  default: assert(false);
  }

  general_assert(samwrite(out, b) >= 0, "Fail to write alignments to BAM file!");
  if (is_paired) general_assert(samwrite(out, b2) >= 0, "Fail to write alignments to BAM file!");

  return true;
}

void BamAlignment::writeToBED(FILE *fo, const bam_header_t *header) {
  if (b->core.flag & 0x0004) return;
  fprintf(fo, "%s\t%d\t%d\t%s\t%.2f\t%c\n", header->target_name[b->core.tid], b->core.pos, b->core.pos + b->core.l_qseq, (char*)bam1_qname(b), getFrac() * 1000.0, ((b->core.flag & 0x0010) == 0 ? '+' : '-'));
  if (!is_paired || (b2->core.flag & 0x0004)) return;
  fprintf(fo, "%s\t%d\t%d\t%s\t%.2f\t%c\n", header->target_name[b2->core.tid], b2->core.pos, b2->core.pos + b2->core.l_qseq, (char*)bam1_qname(b2), getFrac() * 1000.0, ((b2->core.flag & 0x0010) == 0 ? '+' : '-'));
}

void BamAlignment::writeToTagAlign(FILE *fo, const bam_header_t *header) {
  if (b->core.flag & 0x0004) return;
  general_assert(!is_paired, "Paired-end alignments are detected for tagAlign format!");

  SEQstring si;
  assert(getSEQ(si));

  fprintf(fo, "%s\t%d\t%d\t%s\t%.2f\t%c\n", header->target_name[b->core.tid], b->core.pos, b->core.pos + b->core.l_qseq, si.toString().c_str(), getFrac() * 1000.0, ((b->core.flag & 0x0010) == 0 ? '+' : '-'));
}
