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

#ifndef MYHEAP_H_
#define MYHEAP_H_

#include <cassert>
#include <vector>
#include <algorithm>
#include <fstream>
#include <iostream>

#include "utils.h"

struct HeapType {
	int id;
	READ_INT_TYPE nreads;
	HIT_INT_TYPE nlines;
	
	HeapType() { id = -1; nreads = 0; nlines = 0; }
	HeapType(int id, READ_INT_TYPE nreads, HIT_INT_TYPE nlines) { this->id = id; this->nreads = nreads; this->nlines = nlines; }

	HeapType(const HeapType& o) { id = o.id; nreads = o.nreads; nlines = o.nlines; }

	bool operator< (const HeapType& o) const {
		return nlines < o.nlines || (nlines == o.nlines && id < o.id);
	}

	static bool compare(const HeapType& a, const HeapType& b) {
		return a.id < b.id;
	}
};

class MyHeap {
public:
	MyHeap() { size = 0; elements.clear(); }

	void init(int size);
	int getTop() const { return elements[0].id; }
	void updateTop(HIT_INT_TYPE nlines);

	READ_INT_TYPE getNum(int id) const { return elements[id].nreads; }
	void print(const char* outF = NULL);

private:
	int size;
	std::vector<HeapType> elements;
};

inline void MyHeap::init(int size) {
	this->size = size;
	elements.clear();
	for (int i = 0; i < size; i++) 
		elements.push_back(HeapType(i, 0, 0));
}

inline void MyHeap::updateTop(HIT_INT_TYPE nlines) {
	int pos = 0, minpos;
	++elements[pos].nreads;
	elements[pos].nlines += nlines;
	while (pos * 2 + 1 < size) {
		minpos = pos; 
		if (elements[pos * 2 + 1] < elements[pos]) minpos = pos * 2 + 1;
		if ((pos * 2 + 2 < size) && (elements[pos * 2 + 2] < elements[minpos])) minpos = pos * 2 + 2;
		if (pos == minpos) break;
		std::swap(elements[pos], elements[minpos]);
		pos = minpos;
	}
}

// After calling this function, the heap structure cannot be used
inline void MyHeap::print(const char* outF) {
	std::sort(elements.begin(), elements.end(), HeapType::compare);

	if (outF != NULL) {
		std::ofstream fout(outF);
		assert(fout.is_open());
		
		for (int i = 0; i < size; i++) 
			fout<< elements[i].id<< '\t'<< elements[i].nreads<< '\t'<< elements[i].nlines<< std::endl;
	}

	if (verbose) {
		for (int i = 0; i < size; i++) 
			std::cout<< "Partition "<< elements[i].id<< " has "<< elements[i].nreads<< " reads and "<< elements[i].nlines<< " lines."<< std::endl;
	}
}

#endif
