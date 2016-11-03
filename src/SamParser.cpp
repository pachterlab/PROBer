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

#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <cassert>

#include "htslib/sam.h"

#include "my_assert.h"
#include "SamParser.hpp"

SamParser::SamParser(const char* inpF, bam_hdr_t* input_header) {
	sam_in = sam_open(inpF, "r");
	general_assert(sam_in != 0, "Cannot open " + cstrtos(inpF) + "! It may not exist.");

	delete_header = true;
	header = sam_hdr_read(sam_in);
	general_assert(header != 0, "Fail to parse the header!");

	if (input_header != NULL) {
		bam_hdr_destroy(header);
		header = input_header;
		delete_header = false;
	}

	memset(program_id, 0, sizeof(program_id));
}

SamParser::~SamParser() {
	if (delete_header) bam_hdr_destroy(header);
	sam_close(sam_in);
}

// This is an simple implementation, improve it later if necessary
const char* SamParser::getProgramID() {
	if (program_id[0]) return program_id;
	
	char *p = strstr(header->text, "@PG\t");
	assert(p != NULL);
	p += 4;

	char *fr = p;
	while (*p != '\n' && *p != '\0') {
		if (*p == '\t') {
			if (p - fr > 3 && !strncmp(fr, "ID:", 3)) {
	assert(p - fr - 3 <= 100);
	strncpy(program_id, fr + 3, p - fr - 3);
	return program_id;
			}
			fr = p + 1;
		}
		++p;
	}
	if (p - fr > 3 && !strncmp(fr, "ID:", 3)) {
		assert(p - fr - 3 <= 100);
		strncpy(program_id, fr + 3, p - fr - 3);
		return program_id;
	}
	assert(false);
}
