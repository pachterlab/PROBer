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

#ifndef PROBERREADMODEL_H_
#define PROBERREADMODEL_H_

#include <cassert>
#include <string>
#include <fstream>
#include <algorithm>

#include "utils.h"
#include "sampling.hpp"

#include "RefSeq.hpp"
#include "Refs.hpp"

#include "SEQstring.hpp"
#include "QUALstring.hpp"
#include "CIGARstring.hpp"

#include "BamAlignment.hpp"
#include "AlignmentGroup.hpp"

#include "MateLenDist.hpp"
#include "SequencingModel.hpp"
#include "NoiseProfile.hpp"
#include "QualDist.hpp"

#include "InMemoryStructs.hpp"

class PROBerReadModel {
public:

	/*
		@param   model_type   0, SE, no qual; 1, SE qual; 2, PE, no qual; 3 PE, qual
		@param   refs         a pointer to reference sequences
		@param   read_length  if set, assume all read length are the same, in addition, all reads whose lengths < read_length are due to adaptor trimming
		@comment: Master thread for learning model parameters
	 */
	PROBerReadModel(int model_type, Refs* refs, int read_length = -1);

	/*
		@param   master_model   The master model for learning parameters
		@comment: Constructor function for slave threads
	 */
	PROBerReadModel(PROBerReadModel* master_model);

	/*
		@param     refs     a pointer to the reference sequences
		@param     sampler  a pointer to a sampler for simulation
		@comment:  Used for simulation only, read learned parameters from files 
	 */
	PROBerReadModel(Refs* refs, Sampler* sampler);
	
	~PROBerReadModel();

	/*
		@func   get model type
	 */
	int getModelType() const { return model_type; }

	void update_preprocess(AlignmentGroup& ag, bool isAligned);

	void finish_preprocess();

	/*
		@param   ag_in_mem   an in-memory alignment group, recorded information necessary for EM iteration
		@param   aligns      a pointer to all alignments in ag_in_mem
		@param   ag   an alignment group, which contains the read sequence etc.
		@func   set conditional probabilities to an alignments of a read
	 */
	void setConProbs(InMemAlignG* ag_in_mem, InMemAlign* aligns, AlignmentGroup& ag);

	/*
		@param   ag_in_mem
		@param   aligns
		@param   ag
		@param   noise_frac   fractional weight at noise transcript
	 */
	void update(InMemAlignG* ag_in_mem, InMemAlign* aligns, AlignmentGroup& ag, double noise_frac);

	/*
		@return  partial log-likelihood for unalignable reads
	 */
	double calcLogP() {
		return loglik + npro->calcLogP();
	}

	void init();
	void collect(PROBerReadModel* o);
	void finish();

	void read(const char* modelF);
	void write(const char* modelF);

	void simulate(READ_INT_TYPE rid, int tid, int pos, int fragment_length, std::ofstream* out1, std::ofstream* out2 = NULL);

	void startSimulation();
	void finishSimulation();

private:
	int model_type; // 0, SE, no Qual; 1, SE, Qual; 2, PE, no Qual; 3, PE, Qual

	MateLenDist *mld1, *mld2; // mld1, mate length distribution 1; mld2, mate length distribution 2.
	QualDist *qd;
	SequencingModel *seqmodel;
	NoiseProfile *npro;

	int max_len; // maximum mate length
	double loglik; // partial log-likelihood for unaligned reads

	Refs *refs;  
	Sampler *sampler;

	int read_length; // the minimum read length, if read_length is set (not -1), all mates have a same length.
};

inline void PROBerReadModel::update_preprocess(AlignmentGroup& ag, bool isAligned) {
	// Update MLDs
	int len = read_length < 0 ? ag.getSeqLength(1) : read_length;
	mld1->update(len, !isAligned);
	if (model_type >= 2) {
		len = read_length < 0 ? ag.getSeqLength(2) : read_length;
		mld2->update(len, !isAligned);
	}
	
	// Updae QualDist
	if (model_type & 1) {
		QUALstring qual;
		ag.getQUAL(qual); qd->update(qual);
		if (model_type == 3) {
			ag.getQUAL(qual, 2); qd->update(qual);
		}
	}
	
	// Update NoiseProfile
	if (!isAligned) {
		SEQstring seq;
		ag.getSEQ(seq); npro->updateC(seq);
		if (model_type >= 2) {
			ag.getSEQ(seq, 2); npro->updateC(seq);
		}
	}
}

inline void PROBerReadModel::setConProbs(InMemAlignG* ag_in_mem, InMemAlign* aligns, AlignmentGroup& ag) {
	int seqlen;
	SEQstring seq; // seq, qual and cigar must be in each function since we have multiple threads!
	QUALstring qual;
	CIGARstring cigar;
	const RefSeq* refseq = NULL;
	
	// Get read sequences and quality scores
	assert(ag.getSEQ(seq));
	seqlen = read_length < 0 ? ag.getSeqLength() : read_length;
	if (model_type & 1) assert(ag.getQUAL(qual));
	// set noise probability    
	ag_in_mem->noise_conprb = mld1->getProb(seqlen) * npro->getProb(seq);
	// set alignment probabilities
	for (int i = 0; i < ag_in_mem->size; ++i) if (aligns[i].conprb != -1.0) {
		refseq = refs->getRef(aligns[i].tid);
		assert(ag.getAlignment(i)->getCIGAR(cigar));
		aligns[i].conprb = (aligns[i].fragment_length > 0 ? mld1->getProb(seqlen, aligns[i].fragment_length) : mld1->getProb(seqlen)) * \
			seqmodel->getProb('+', aligns[i].pos, refseq, &cigar, &seq, ((model_type & 1) ? &qual : NULL));
	}
	
	if (model_type >= 2) {
		// paired-end reads
		assert(ag.getSEQ(seq, 2));
		seqlen = read_length < 0 ? ag.getSeqLength(2) : read_length;
		if (model_type & 1) assert(ag.getQUAL(qual, 2));
		ag_in_mem->noise_conprb *= mld2->getProb(seqlen) * npro->getProb(seq);
		for (int i = 0; i < ag_in_mem->size; ++i) if (aligns[i].conprb != -1.0) {
			refseq = refs->getRef(aligns[i].tid);
			assert(ag.getAlignment(i)->getCIGAR(cigar, 2));
			assert(aligns[i].fragment_length > 0);
			aligns[i].conprb *= mld2->getProb(seqlen, aligns[i].fragment_length) * \
	seqmodel->getProb('-', refseq->getLen() - aligns[i].pos - aligns[i].fragment_length, refseq, &cigar, &seq, ((model_type & 1) ? &qual : NULL));
		}
	}
}

inline void PROBerReadModel::update(InMemAlignG* ag_in_mem, InMemAlign* aligns, AlignmentGroup& ag, double noise_frac) {
	SEQstring seq;
	QUALstring qual;
	CIGARstring cigar;
	const RefSeq* refseq = NULL;
	
	assert(ag.getSEQ(seq));
	if (model_type & 1) assert(ag.getQUAL(qual));
	// update noise prob
	npro->update(seq, noise_frac);
	// update alignment probs
	for (int i = 0; i < ag_in_mem->size; ++i) if (aligns[i].frac > 0.0) {
		refseq = refs->getRef(aligns[i].tid);
		assert(ag.getAlignment(i)->getCIGAR(cigar));
		seqmodel->update(aligns[i].frac, '+', aligns[i].pos, refseq, &cigar, &seq, ((model_type & 1) ? &qual : NULL));
	}

	if (model_type >= 2) {
		// paired-end reads
		assert(ag.getSEQ(seq, 2));
		if (model_type & 1) assert(ag.getQUAL(qual, 2));
		// update noise prob
		npro->update(seq, noise_frac);
		// update alignment probs
		for (int i = 0; i < ag_in_mem->size; ++i) if (aligns[i].frac > 0.0) {
			refseq = refs->getRef(aligns[i].tid);
			assert(ag.getAlignment(i)->getCIGAR(cigar, 2));
			assert(aligns[i].fragment_length > 0);
			seqmodel->update(aligns[i].frac, '-', refseq->getLen() - aligns[i].pos - aligns[i].fragment_length, refseq, &cigar, &seq, ((model_type & 1) ? &qual : NULL));
		}
	}
} 

inline void PROBerReadModel::simulate(READ_INT_TYPE rid, int tid, int pos, int fragment_length, std::ofstream* out1, std::ofstream* out2) {
	int m2pos;
	int mateL1, mateL2;
	std::string qual1, qual2, cigar1, cigar2, readseq1, readseq2;
	
	// Simulate reads
	if (tid == 0) {
		cigar1 = cigar2 = "*";

		mateL1 = mld1->simulate(sampler);
		if (model_type & 1) qd->simulate(sampler, mateL1, qual1);
		npro->simulate(sampler, mateL1, readseq1);
		
		if (model_type >= 2) {
			mateL2 = mld2->simulate(sampler);
			if (model_type & 1) qd->simulate(sampler, mateL2, qual2);
			npro->simulate(sampler, mateL2, readseq2);
		}
	}
	else {
		const RefSeq* ref = refs->getRef(tid);
		mateL1 = mld1->simulate(sampler, fragment_length);
		if (model_type & 1) qd->simulate(sampler, mateL1, qual1);
		seqmodel->simulate(sampler, mateL1, '+', pos, ref, qual1, cigar1, readseq1);

		if (model_type >= 2) {
			mateL2 = mld2->simulate(sampler, fragment_length); 
			if (model_type & 1) qd->simulate(sampler, mateL2, qual2);
			m2pos = ref->getLen() - pos - fragment_length;
			seqmodel->simulate(sampler, mateL2, '-', m2pos, ref, qual2, cigar2, readseq2);
		}
	}

	// Output reads
	(*out1)<< ((model_type & 1) ? '@' : '>') << rid<< '_'<< tid<< '_'<< pos<< '_'<< fragment_length<< '_'<< cigar1<< (model_type >= 2 ? "/1" : "") << std::endl;
	(*out1)<< readseq1<< std::endl;
	if (model_type & 1) (*out1)<< '+'<< std::endl<< qual1<< std::endl;

	if (model_type >= 2) {
		(*out2)<< ((model_type & 1) ? '@' : '>') << rid<< '_'<< tid<< '_'<< pos<< '_'<< fragment_length<< '_'<< cigar2<< "/2"<< std::endl;
		(*out2)<< readseq2<< std::endl;
		if (model_type & 1) (*out2)<< '+'<< std::endl<< qual2<< std::endl;
	}
}

#endif
