/*
 * transcripts are numbered from 1. 0 is reserved for noise isoform
 */
#ifndef TRANSCRIPTS_H_
#define TRANSCRIPTS_H_

#include<cassert>
#include<string>
#include<vector>
#include<algorithm>

#include "Transcript.hpp"

class Transcripts {
public:
	Transcripts(int type = 0) {
		M = 0; this->type = type;
		transcripts.clear();
		transcripts.push_back(Transcript());

		e2i = i2e = NULL;
	}

  ~Transcripts() { 
    if (e2i != NULL) delete[] e2i;
    if (i2e != NULL) delete[] i2e;
  }

	int getM() const { return M; }

	// used in shrinking the transcripts
	void setM(int M) { this->M = M; transcripts.resize(M + 1); } 
	
	void move(int from, int to) {
	  assert(from >= to);
	  if (from > to) transcripts[to] = transcripts[from];
	}
	
	int getType() const { return type; }
	void setType(int type) { this->type = type; }

	bool isAlleleSpecific() const { return type == 2; }

	const Transcript& getTranscriptAt(int pos) const {
		assert(pos > 0 && pos <= M);
		return transcripts[pos];
	}

	void add(const Transcript& transcript) {
		transcripts.push_back(transcript);
		++M;
	}

	void sort() {
		std::sort(transcripts.begin(), transcripts.end());
	}

	void readFrom(const char*);
	void writeTo(const char*);

	//Eid: external sid, 0 <= eid < M
	int getInternalSid(int eid) const {
		assert(eid >= 0 && eid < M);
		return e2i[eid];
	}

	// Iid: internal sid, 0 < iid <= M
	int getExternalSid(int iid) const {
	  assert(iid > 0 && iid <= M);
	  return i2e[iid];
	}

  const int* getE2IArray() const {
    return e2i;
  }

	const Transcript& getTranscriptViaEid(int eid) const {
		return transcripts[getInternalSid(eid)];
	}

	void buildMappings(const char* imdName, int n_targets = 0, char** target_name = NULL);

private:
	int M, type; // type 0 from genome, 1 standalone transcriptome, 2 allele-specific 
	std::vector<Transcript> transcripts;

	int *e2i, *i2e; // external sid to internal sid, internal sid to external sid
};

#endif /* TRANSCRIPTS_H_ */
