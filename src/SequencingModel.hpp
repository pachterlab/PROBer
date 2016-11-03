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

#ifndef SEQUENCINGMODEL_H_
#define SEQUENCINGMODEL_H_

#include <cmath>
#include <cassert>
#include <fstream>

#include "RefSeq.hpp"
#include "CIGARstring.hpp"
#include "SEQstring.hpp"
#include "QUALstring.hpp"

#include "Markov.hpp"
#include "Profile.hpp"
#include "QProfile.hpp"

/**
 * We have two explanation of indel models. One is to assume indels are due to difference in the reference sequence and real sequence. The other assumes that sequencing errors cause indels.
 *
 * Using the first assumption, indels are generated from reference and point errors are generated from sequencing process. However, this assumption makes it hard to determine insert size.
 * In addition, different reads generated from a same position may represent different reference sequences, which is also confusing.
 * 
 * With the second assumption, we need to argue why profile and qprofile can still use their position and quality scores for estimation with the present of indels. However, it is consistent 
 * with insert size. Thus we choose the second interpretation.
 */

class SequencingModel {
public:
	SequencingModel(bool hasQual, int maxL = 1000);
	~SequencingModel();

	double getProb(char dir, int pos, const RefSeq* refseq, const CIGARstring* cigar, const SEQstring* seq, const QUALstring* qual = NULL);
	void update(double frac, char dir, int pos, const RefSeq* refseq, const CIGARstring* cigar, const SEQstring* seq, const QUALstring* qual = NULL);

	void init();
	void collect(const SequencingModel* o);
	void finish();

	void read(std::ifstream& fin);
	void write(std::ofstream& fout);

	void simulate(Sampler *sampler, int len, char dir, int pos, const RefSeq* refseq, const std::string& qual, std::string& cigar, std::string& seq);

	void startSimulation();
	void finishSimulation();

private:
	bool hasQual; // if quality scores are available

	Markov *markov;
	Profile *profile;
	QProfile *qprofile;

	void push(std::string& cigar, char opchr, int oplen) {
		int s = 0;
		char arr[50];

		while (oplen > 0) { arr[s++] = oplen % 10 + '0'; oplen /= 10; }
		while (s > 0) cigar.push_back(arr[--s]);
		cigar.push_back(opchr);
	}
};

/*
	@param   dir      '+' or '-'  
	@param   pos      position in its own strand, 0-based
	@param   refseq   reference sequence in '+' strand
	@param   cigar    cigar string
	@param   seq      read sequence
	@param   qual     quality score string
	@return  probability of generating such a read
 */
inline double SequencingModel::getProb(char dir, int pos, const RefSeq* refseq, const CIGARstring* cigar, const SEQstring* seq, const QUALstring* qual) {
	double prob = 1.0;
	int len = cigar->getLen();
	int readpos = 0;

	char opchr, last_opchr = 0;
	int oplen;

	int ref_base, read_base;

	for (int i = 0; i < len; ++i) {
		opchr = cigar->opchrAt(i);
		oplen = cigar->oplenAt(i);

		// Markov model probabilities
		prob *= (last_opchr == 0 ? markov->getProb(opchr) : markov->getProb(last_opchr, opchr));
		if (oplen > 1) prob *= pow(markov->getProb(opchr, opchr), oplen - 1);

		if (opchr == 'M' || opchr == '=' || opchr == 'X') {
			for (int j = 0; j < oplen; ++j) {
	ref_base = refseq->baseCodeAt(dir, pos);
	read_base = seq->baseCodeAt(readpos);
	prob *= (hasQual ? qprofile->getProb(qual->qualAt(readpos), ref_base, read_base) : profile->getProb(readpos, ref_base, read_base));
	++pos; ++readpos;
			}
		}
		else if (opchr == 'I') {
			for (int j = 0; j < oplen; ++j) {
	prob *= markov->getIBaseProb(seq->baseCodeAt(readpos));
	++readpos;
			}
		}
		else {
			assert(opchr == 'D');
			pos += oplen;
		}

		last_opchr = opchr;
	}

	return prob;
}

/*
	@param   frac     fractional weight
	@param   dir      '+' or '-'  
	@param   pos      position in its own strand, 0-based
	@param   refseq   reference sequence in '+' strand
	@param   cigar    cigar string
	@param   seq      read sequence
	@param   qual     quality score string
 */
inline void SequencingModel::update(double frac, char dir, int pos, const RefSeq* refseq, const CIGARstring* cigar, const SEQstring* seq, const QUALstring* qual) {
	int len = cigar->getLen();
	int readpos = 0;

	char opchr, last_opchr = 0;
	int oplen;

	int ref_base, read_base;

	for (int i = 0; i < len; ++i) {
		opchr = cigar->opchrAt(i);
		oplen = cigar->oplenAt(i);

		// Markov model probabilities
		if (last_opchr == 0) markov->update(opchr, frac);
		else markov->update(last_opchr, opchr, frac);

		if (oplen > 1) markov->update(opchr, opchr, frac * (oplen - 1));

		if (opchr == 'M' || opchr == '=' || opchr == 'X') {
			for (int j = 0; j < oplen; ++j) {
	ref_base = refseq->baseCodeAt(dir, pos);
	read_base = seq->baseCodeAt(readpos);
	if (hasQual) qprofile->update(qual->qualAt(readpos), ref_base, read_base, frac);
	else profile->update(readpos, ref_base, read_base, frac);
	++pos; ++readpos;
			}
		}
		else if (opchr == 'I') {
			for (int j = 0; j < oplen; ++j) {
	markov->updateIBase(seq->baseCodeAt(readpos), frac);
	++readpos;
			}
		}
		else {
			assert(opchr == 'D');
			pos += oplen;
		}

		last_opchr = opchr;
	}
}

/*
	@param   sampler   Sampler used for simulation
	@param   len       length of the read
	@param   dir       strand of the read, either '+' or '-'
	@param   pos
	@param   refseq
	@param   qual      a string of (SANGER score + 33)
	@param   cigar
	@param   seq
	@comment: If the end of sequence is reached, only insertions are generated
 */
inline void SequencingModel::simulate(Sampler *sampler, int len, char dir, int pos, const RefSeq* refseq, const std::string& qual, std::string& cigar, std::string& seq) {
	int readpos = 0;
	char opchr, last_opchr;
	int oplen;
	int refLen = refseq->getLen();

	cigar.clear();
	seq.clear();

	last_opchr = 0; oplen = 0;
	while (readpos < len) {
		if (pos >= refLen) opchr = 'I';
		else if (last_opchr == 0) opchr = markov->simulate(sampler);
		else opchr = markov->simulate(sampler, last_opchr);
		
		if (opchr == 'M') {
			seq.push_back(hasQual ? qprofile->simulate(sampler, qual[readpos] - 33, refseq->baseCodeAt(dir, pos)) : profile->simulate(sampler, readpos, refseq->baseCodeAt(dir, pos)));
			++readpos; ++pos;
		}
		else if (opchr == 'I') {
			seq.push_back(markov->simulateIBase(sampler));
			++readpos;
		}
		else {
			assert(opchr == 'D');
			++pos;
		}

		if (last_opchr != opchr) { 
			if (last_opchr != 0) push(cigar, last_opchr, oplen);
			last_opchr = opchr;
			oplen = 0;
		}
		++oplen;
	}
	if (last_opchr != 0) push(cigar, last_opchr, oplen);
}

#endif
