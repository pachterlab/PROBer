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

#ifndef BAMWRITER_H_
#define BAMWRITER_H_

#include <string>

#include "htslib/sam.h"

#include "BamAlignment.hpp"
#include "AlignmentGroup.hpp"

class BamWriter {
public:
	// if program_id is NULL, use the full header
	BamWriter(const char* outF, const bam_hdr_t* header, const char* program_id = NULL);
	~BamWriter();

	// choice: 0, do nothing; 1, delete read sequence and qual score; 2, add read sequence and qual score  
	bool write(BamAlignment& b, int choice = 0) { return b.write(bam_out, header, choice); }
	
	// choice: 0, do nothing; 1, delete read sequence and qual score; 2, add read sequence and qual score
	bool write(AlignmentGroup& ag, int choice = 0) { return ag.write(bam_out, header, choice); } 

private:
	samFile* bam_out;
	bam_hdr_t* header;
	
	bam_hdr_t* header_duplicate_without_text(const bam_hdr_t* ori_h);
	void header_append_new_text(bam_hdr_t* header, const std::string& new_text);
};

#endif
