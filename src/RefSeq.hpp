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

#ifndef REFSEQ_H_
#define REFSEQ_H_

#include <cctype>
#include <cassert>
#include <fstream>
#include <string>

#include <stdint.h>

#include "utils.h"
#include "my_assert.h"
#include "CIGARstring.hpp"
#include "MDstring.hpp"
#include "SEQstring.hpp"

class RefSeq {
public:
	RefSeq();
	RefSeq(const std::string& name, const std::string& rawseq);
	RefSeq(const RefSeq& o);
	RefSeq& operator= (const RefSeq&);
	
	bool read(std::ifstream& fin);
	void write(std::ofstream& fout);

	int getLen() const { return len; }
	const std::string& getName() const { return name; }
	const std::string& getSeq() const { return seq; }

	char baseAt(char dir, int pos) const {
		assert(pos >= 0 && pos < len);
		return (dir == '+' ? seq[pos] : base2rbase[seq[len - pos - 1]]);
	}

	int baseCodeAt(char dir, int pos) const {
		assert(pos >= 0 && pos < len);
		return (dir == '+' ? base2code[seq[pos]] : rbase2code[seq[len - pos - 1]]);
	}

	void setUp(char dir, CIGARstring& cigar, MDstring& mdstr, SEQstring& seq);
	
private:
	int len; // the transcript length
	std::string name; // the tag
	std::string seq; // the sequence, in the forward strand

	void convertRawSeq() {
		for (int i = 0; i < len; ++i) {
			general_assert(isalpha(seq[i]), "Sequence contains unknown code " + ctos(seq[i]) + "!");
			seq[i] = toupper(seq[i]);
			if (seq[i] != 'A' && seq[i] != 'C' && seq[i] != 'G' && seq[i] != 'T') seq[i] = 'N';
		}
	}
	
};

/*
	@function constructing RefSeq based on MD filed
	@param   dir     alignment direction
	@param   cigar   CIGAR string
	@param   mdstr   MD string
	@param   seqstr  SEQ string
 */
inline void RefSeq::setUp(char dir, CIGARstring& cigar, MDstring& mdstr, SEQstring& seqstr) {
	len = 0; name = seq = "";

	char old_dir = cigar.getDir();
	cigar.setDir(dir);
	seqstr.setDir(dir);
	
	int pos = -1, cigar_len = cigar.getLen();
	int optype, oplen;
	char ref_base;
	
	for (int i = 0; i < cigar_len; ++i) {
		optype = cigar.optypeAt(i);
		oplen = cigar.oplenAt(i);
		for (int j = 0; j < oplen; ++j) {
			if (optype & 1) ++pos;
			if (optype & 2) {
	ref_base = mdstr.next();
	assert(ref_base >= 0);
	if (ref_base == 0) ref_base = seqstr.baseAt(pos);
	seq += ref_base;
			}
		}
	}
	
	name = "constructed_from_BAM_alignment";
	len = seq.length();

	cigar.setDir(old_dir);
	seqstr.setDir(old_dir);
}

#endif
