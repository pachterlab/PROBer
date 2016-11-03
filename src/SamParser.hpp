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

#ifndef SAMPARSER_H_
#define SAMPARSER_H_

#include <string>

#include "htslib/sam.h"

#include "BamAlignment.hpp"
#include "AlignmentGroup.hpp"

class SamParser {
public:
	SamParser(const char* inpF, bam_hdr_t* input_header = NULL);
	~SamParser();

	const bam_hdr_t* getHeader() const { return header; }

	bam_hdr_t* pass_header() { delete_header = false; return header; } // pass the header to an outside variable, then it is the outside variable's responsibility to delete header

	const char* getProgramID(); // scan header to look up program ID, slow

	bool next(BamAlignment& b) { return b.read(sam_in, header); }

	bool next(AlignmentGroup& ag) { return ag.read(sam_in, header); }
	
private:
	samFile* sam_in;
	bam_hdr_t* header;

	char program_id[1005];
	bool delete_header;
};

#endif

