/* Copyright (c) 2015
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

#include<cassert>
#include<string>
#include<fstream>

#include "CIGARstring.hpp"
#include "MDstring.hpp"
#include "SEQstring.hpp"
#include "RefSeq.hpp"

RefSeq::RefSeq() {
  len = 0;
  name = seq = "";
}

/*
  @function   constructor.
  @param   name   transcript name
  @param   rawseq    raw transcript sequence
 */
RefSeq::RefSeq(const std::string& name, const std::string& rawseq) {
  len = rawseq.length();
  assert(len > 0);
  
  this->name = name;
  this->seq = rawseq;
  
  convertRawSeq();
}

RefSeq::RefSeq(const RefSeq& o) {
  len = o.len;
  name = o.name;
  seq = o.seq;
}

RefSeq& RefSeq::operator= (const RefSeq &rhs) {
  if (this != &rhs) {
    len = rhs.len;
    name = rhs.name;
    seq = rhs.seq;
  }
  
  return *this;
}

bool RefSeq::read(std::ifstream& fin) {
  std::string line;

  if (!getline(fin, name)) return false;
  name = name.substr(1);
  if (!getline(fin, seq)) return false;
  len = seq.length();

  return true;
}

void RefSeq::write(std::ofstream& fout) {
  fout<< ">"<< name<< std::endl;
  fout<< seq<< std::endl;
}

/*
  @function constructing RefSeq based on MD filed
  @param   dir     alignment direction
  @param   cigar   CIGAR string
  @param   mdstr   MD string
  @param   seq     SEQ string
 */
void RefSeq::setUp(char dir, CIGARstring& cigar, MDstring& mdstr, SEQstring& seq) {
  len = 0; name = seq = "";

  char old_dir = cigar.getDir();
  cigar.setDir(dir);
  seq.setDir(dir);
  
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
	if (ref_base == 0) ref_base = seq.baseAt(pos);
	seq += ref_base;
      }
    }
  }
  
  name = "constructed_from_BAM_alignment";
  len = seq.length();

  cigar.setDir(old_dir);
  seq.setDir(old_dir);
}
