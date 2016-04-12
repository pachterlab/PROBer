#include<cmath>
#include<cstdio>
#include<cstring>
#include<cassert>
#include<string>
#include<algorithm>

#include<stdint.h>
#include "htslib/sam.h"

#include "my_assert.h"
#include "BamAlignment.hpp"


const uint8_t BamAlignment::rnt_table[16] = {0, 8, 4, 0, 2, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 15};

bool BamAlignment::read(samFile* in, bam_hdr_t* header, BamAlignment* o) {
  is_aligned = -1;

  if (b == NULL) b = bam_init1();
  if (sam_read1(in, header, b) < 0) return false;

  is_paired = bam_is_paired(b);
  // read the second mate
  if (is_paired) { 
    if (b2 == NULL) b2 = bam_init1();

    general_assert(sam_read1(in, header, b2) >= 0 && bam_is_paired(b2), "Fail to read the other mate for a paired-end alignment!");
    general_assert(((b->core.flag & 0x00C0) == 0x0040 && (b2->core.flag & 0x00C0) == 0x0080) || 
		   ((b->core.flag & 0x00C0) == 0x0080 && (b2->core.flag & 0x00C0) == 0x0040), 
		   "Cannot detect both mates of a paired-end alignment!");
    
    if (bam_is_read2(b)) { bam1_t *tmp = b; b = b2; b2 = tmp; }
  }

  // calculate is_aligned
  is_aligned = bam_is_mapped(b);
  if (is_paired) is_aligned |= ((char)bam_is_mapped(b2) << 1);
  
  // The following four statements are grouped together
  if (b->core.l_qseq <= 0) b->core.l_qseq = o->getSeqLength(1);
  if (is_paired && b2->core.l_qseq <= 0) b2->core.l_qseq = o->getSeqLength(2);
  assert(!(is_aligned & 1) || b->core.l_qseq == bam_cigar2qlen(b->core.n_cigar, bam_get_cigar(b)));
  assert(!(is_aligned & 2) || b2->core.l_qseq == bam_cigar2qlen(b2->core.n_cigar, bam_get_cigar(b2)));

  return true;
}

// choice: 0, do nothing; 1, delete read sequence and qual score; 2, add read sequence and qual score. o, the alignment that contain read sequence and quality score information
bool BamAlignment::write(samFile* out, const bam_hdr_t* header, int choice, BamAlignment* o) {
  assert(is_aligned >= 0 && b != NULL && (!is_paired || b2 != NULL));

  if (b->core.l_qname == 1) b->core.l_qseq = 0;
  if (is_paired && (b2->core.l_qname == 1)) b2->core.l_qseq = 0;

  switch(choice) {
  case 0: 
    break;
  case 1:
    if (b->core.l_qname > 1) compress(b);
    if (is_paired && (b2->core.l_qname > 1)) compress(b2);
    break;
  case 2:
    if (b->core.l_qname == 1) decompress(b, o->b);
    if (is_paired && (b2->core.l_qname == 1)) decompress(b2, o->b2);
    break;
  default: assert(false);
  }

  general_assert(sam_write1(out, header, b) >= 0, "Fail to write alignments to BAM file!");
  if (is_paired) general_assert(sam_write1(out, header, b2) >= 0, "Fail to write alignments to BAM file!");

  return true;
}
