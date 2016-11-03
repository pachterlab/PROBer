/* Copyright (c) 2016
   Bo Li (University of California, Berkeley)
   bli25@berkeley.edu

   This program is free software; you can redistribute it and/or
   modify it under the terms of the GNU General Public License as
   published by the Free Software Foundation; either version 3 of the
   License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
   General Public License for more details.   

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
   USA
*/

#ifndef BAMALIGNMENT_H_
#define BAMALIGNMENT_H_

#include <cmath>
#include <cassert>
#include <string>
#include <algorithm>

#include <stdint.h>
#include "htslib/sam.h"

#include "my_assert.h"

#include "CIGARstring.hpp"
#include "SEQstring.hpp"
#include "QUALstring.hpp"
#include "MDstring.hpp"

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
	
	/*
		@param   in       input SAM/BAM/CRAM file handler
		@param   header   input file header // consider change bam_hdr_t* to const bam_hdr_t* later
		@param   o        optional BAM alignment
	 */
	bool read(samFile* in, bam_hdr_t* header, BamAlignment* o = NULL);

	/*
		@param   out      output BAM file handler
		@param   header   output file header
		@param   choice   0, do nothing; 1, delete read sequence and qual score; 2, add read sequence and qual score
		@param   o        optional BAM alignment
	 */
	bool write(samFile* out, const bam_hdr_t* header, int choice = 0, BamAlignment* o = NULL);
	
	// overall stats
	
	bool isPaired() const { return is_paired; }
	
	// -1: nothing is loaded, 0: not aligned, 1: first mate is aligned (including SE reads), 2: second mate is aligned, 3: fully aligned 
	int isAligned() const { return is_aligned; }

	int getInsertSize() const {
		assert(is_aligned == 3);
		return abs(b->core.isize); 
	}
	
	int getTid() const { 
		assert(is_aligned > 0); 
		return (is_aligned & 1) ? b->core.tid : b2->core.tid;
	}
	
	char getDir() const { 
		assert(is_aligned > 0);
		if (is_aligned & 1) return !bam_is_rev(b) ? '+' : '-';
		return bam_is_rev(b2) ? '+' : '-';
	}

	char getMateDir(int mate = 1) const {
		assert((mate == 1 && bam_is_mapped(b)) || (mate == 2 && is_paired && bam_is_mapped(b2)));
		if (mate == 1) return bam_is_rev(b) ? '-' : '+';
		else return bam_is_rev(b2) ? '-' : '+';
	}

	// designed for iCLIP (mate == 1) / eCLIP (mate == 2), return cross link site in forward strand coordinate
	int getCrosslinkSite(int mate = 1) const {
		if (mate == 1) return bam_is_rev(b) ? bam_endpos(b) : b->core.pos - 1;
		else return bam_is_rev(b2) ? bam_endpos(b2) : b2->core.pos - 1;
	}
	
	/*
		@param     fragment_length     The average fragment length, 0 means no fragment length is provided
		@return    If fragment_length == 0, return the leftmost position of two mates. Otherwise, return the leftmost position calculated with fragment length.
							 If the calculated position < 0, set it to 0
		@comment   bam_endpos gives the right-most position + 1 
	 */
	int getLeftMostPos(int fragment_length = 0) const {
		assert(is_aligned > 0);
		if (is_aligned == 3) return std::min(b->core.pos, b2->core.pos);
		if (fragment_length <= 0) return (is_aligned & 1) ? b->core.pos : b2->core.pos;

		if (is_aligned & 1) return bam_is_rev(b) ? std::max(0, bam_endpos(b) - fragment_length) : b->core.pos;
		else return bam_is_rev(b2) ? std::max(0, bam_endpos(b2) - fragment_length) : b2->core.pos;
	}
	
	// if length = 2k, midpoint is k - 1
	int getMidPos(int fragment_length = 0) const {
		assert(is_aligned > 0);
		return (is_aligned == 3 ? getLeftMostPos() + (getInsertSize() - 1) / 2 : getLeftMostPos(fragment_length) + (fragment_length - 1) / 2);
	}
	
	int getMapQ() const { return (is_aligned & 1) ? b->core.qual : b2->core.qual; }
	
	void setMapQ(int MapQ) { 
		b->core.qual = (bam_is_mapped(b) ? MapQ : 255);
		if (is_paired) b2->core.qual = (bam_is_mapped(b2) ? MapQ : 255);
	}
	
	// for mates

	/*
		@param   mate   which mate (1 or 2)
		@param   target_len   the length of reference sequence
		@comment: This function returns the smallest position of the mate from its strand
	 */
	int getDirPos(int mate, int target_len) const {
		assert(is_aligned > 0 && (mate == 1 || (is_paired && mate == 2)));
		if (mate == 1) return bam_is_rev(b) ? target_len - bam_endpos(b) : b->core.pos;
		else return bam_is_rev(b2) ? target_len - bam_endpos(b2) : b2->core.pos;
	}

	// mate = 0, the name of the read; 1, mate 1; 2, mate 2
	const char* getName(int mate = 0) const { 
		return bam_get_qname(mate < 2 ? b : b2); 
	}
	
	// length of the query sequence
	int getSeqLength(int mate = 1) const { 
		assert(mate == 1 || (is_paired && mate == 2));
		assert(((mate == 1) && (b->core.l_qseq > 0)) || ((mate == 2) && (b2->core.l_qseq > 0)));
		return mate == 1 ? b->core.l_qseq : b2->core.l_qseq;
	} 
	
	int getAlignedLength(int mate = 1) const {
		assert(mate == 1 || (is_paired && mate == 2));
		return mate == 1 ? bam_cigar2rlen(b->core.n_cigar, bam_get_cigar(b)) : bam_cigar2rlen(b2->core.n_cigar, bam_get_cigar(b2));
	}

	bool getCIGAR(CIGARstring& ci, int mate = 1) {
		assert(mate == 1 || (is_paired && mate == 2));
		if (mate == 1) ci.setUp(bam_get_cigar(b), b->core.n_cigar, is_ori(b));
		else ci.setUp(bam_get_cigar(b2), b2->core.n_cigar, is_ori(b2));
		return true;
	}
	
	bool getSEQ(SEQstring& si, int mate = 1) {
		assert(mate == 1 || (is_paired && mate == 2));
		assert(((mate == 1) && (b->core.l_qseq > 0)) || ((mate == 2) && (b2->core.l_qseq > 0)));
		if (mate == 1) si.setUp(bam_get_seq(b), b->core.l_qseq, is_ori(b));
		else si.setUp(bam_get_seq(b2), b2->core.l_qseq, is_ori(b2));
		return true;
	}
	
	bool getQUAL(QUALstring& qi, int mate = 1) {
		assert(mate == 1 || (is_paired && mate == 2));
		assert(((mate == 1) && (b->core.l_qseq > 0)) || ((mate == 2) && (b2->core.l_qseq > 0)));
		if (mate == 1) qi.setUp(bam_get_qual(b), b->core.l_qseq, is_ori(b));
		else qi.setUp(bam_get_qual(b2), b2->core.l_qseq, is_ori(b2));
		return true;
	}
	
	bool getMD(MDstring& mdstr, int mate = 1) {
		uint8_t* p;    
		assert(mate == 1 || (is_paired && mate == 2));
		p = bam_aux_get((mate == 1 ? b : b2), "MD");
		assert(p != NULL && *p++ == 'Z');
		mdstr.setUp((char*)p);
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
		@func   return true if ZF tag is detected!
	 */
	bool isFiltered() const {
		return bam_aux_get(b, "ZF") != NULL;
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
	double tag2d(uint8_t* p) { return bam_aux2f(p); }
	char* tag2Z(uint8_t* p) { return bam_aux2Z(p); }
	char* tag2H(uint8_t* p) { return bam_aux2Z(p); }
	char* tag2B(uint8_t* p) { return (char*)(p + 1); }

	// If no ZW field, return -1.0
	float getFrac() {
		uint8_t *p_tag = bam_aux_get(((is_aligned & 1) ? b : b2), "ZW");
		return p_tag != NULL ? *(float*)(p_tag + 1) : -1.0;
	}
	
	void setFrac(float frac) {
		uint8_t *p_tag = NULL;

		if (bam_is_mapped(b)) {
			p_tag = bam_aux_get(b, "ZW");
			if (p_tag != NULL) {
	memcpy(p_tag + 1, (uint8_t*)&(frac), bam_aux_type2size('f'));
			} else {
	bam_aux_append(b, "ZW", 'f', bam_aux_type2size('f'), (uint8_t*)&frac);
			}
		}
		
		if (is_paired && bam_is_mapped(b2)) {
			p_tag = bam_aux_get(b2, "ZW");
			if (p_tag != NULL) {
	memcpy(p_tag + 1, (uint8_t*)&(frac), bam_aux_type2size('f'));
			} else {
	bam_aux_append(b2, "ZW", 'f', bam_aux_type2size('f'), (uint8_t*)&frac);
			}
		}

		setMapQ(frac2MapQ(frac));
	}

protected:
	static const uint8_t rnt_table[16];
	
	bool is_paired;
	char is_aligned; // 2 bits, from right to left, the first bit represents first mate and the second bit represents the second mate
									 // Thus, 0, unalignable; 1, only first mate; 2, only second mate; 3, both mates
	bam1_t *b, *b2;


	bool bam_is_paired(const bam1_t* b) const { return (b->core.flag & BAM_FPAIRED); }
	bool bam_is_proper(const bam1_t* b) const { return (b->core.flag & BAM_FPROPER_PAIR); }
	bool bam_is_mapped(const bam1_t* b) const { return !(b->core.flag & BAM_FUNMAP); }
	bool bam_is_unmapped(const bam1_t* b) const { return (b->core.flag & BAM_FUNMAP); }
	bool bam_is_read1(const bam1_t* b) const { return (b->core.flag & BAM_FREAD1); }
	bool bam_is_read2(const bam1_t* b) const { return (b->core.flag & BAM_FREAD2); }
	
	int bam_aux_type2size(char x) {
		if (x == 'C' || x == 'c' || x == 'A') return 1;
		else if (x == 'S' || x == 's') return 2;
		else if (x == 'I' || x == 'i' || x == 'f') return 4;
		else if (x == 'd') return 8;
		else return 0;
	}
	
	uint8_t frac2MapQ(float val) {
		float err = 1.0 - val;
		if (err <= 1e-10) return 100;
		return (uint8_t)(-10 * log10(err) + .5); // round it
	}

	bool is_ori(bam1_t* b) {
		return ((bam_is_unmapped(b) || !bam_is_rev(b)) ? true : false);
	}

	void compress(bam1_t* b) {
		int l_aux = bam_get_l_aux(b);
		b->data[0] = 0;
		memmove(b->data + 1, bam_get_cigar(b), b->core.n_cigar * 4);
		memmove(b->data + 1 + b->core.n_cigar * 4, bam_get_aux(b), l_aux);
		b->l_data = 1 + b->core.n_cigar * 4 + l_aux;
		b->core.l_qname = 1;
		b->core.l_qseq = 0;
	}
	
	// Caution: this function may change b->data's adddress!  
	void expand_data_size(bam1_t* b) {
		if (b->m_data < b->l_data) {
			b->m_data = b->l_data;
			kroundup32(b->m_data);
			b->data = (uint8_t*)realloc(b->data, b->m_data);
		}
	}

	void copy_rc_seq(uint8_t* dst, uint8_t* src, int len) {
		uint8_t base;
		for (int i = 0; i < len; ++i) {
			base = rnt_table[bam_seqi(src, len - i - 1)];
			assert(base > 0);
			if (i & 1) { *dst |= base; ++dst; }
			else { *dst = base << 4; }
		}
	}

	void copy_r_qual(uint8_t* dst, uint8_t* src, int len) {
		for (int i = 0; i < len; ++i) dst[i] = src[len - i - 1];
	}
	
	void decompress(bam1_t* b, bam1_t* other) {
		int l_aux = bam_get_l_aux(b);
		b->core.l_qname = other->core.l_qname;
		b->core.l_qseq = other->core.l_qseq;
		b->l_data = b->core.l_qname + b->core.n_cigar * 4 + (b->core.l_qseq + 1) / 2 + b->core.l_qseq + l_aux;
		expand_data_size(b);
		memmove(bam_get_aux(b), b->data + 1 + b->core.n_cigar * 4, l_aux); // move aux options
		memmove(bam_get_cigar(b), b->data + 1, b->core.n_cigar * 4); // move cigar string
		memcpy(bam_get_qname(b), bam_get_qname(other), b->core.l_qname); // copy qname

		if (bam_is_rev(b) == bam_is_rev(other)) {
			memcpy(bam_get_seq(b), bam_get_seq(other), (b->core.l_qseq + 1) / 2); // copy seq
			memcpy(bam_get_qual(b), bam_get_qual(other), b->core.l_qseq); // copy qual
		}
		else {
			copy_rc_seq(bam_get_seq(b), bam_get_seq(other), b->core.l_qseq); // copy reverse complement seq
			copy_r_qual(bam_get_qual(b), bam_get_qual(other), b->core.l_qseq); // copy reverse qual
		}
	}
};

#endif
