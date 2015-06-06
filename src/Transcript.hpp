#ifndef TRANSCRIPT_H_
#define TRANSCRIPT_H_

#include<string>
#include<vector>
#include<fstream>

/**
   If no genome is provided, seqname field is used to store the allele name.
 */

struct Interval {
	int start, end;

	Interval(int start, int end) {
		this->start = start;
		this->end = end;
	}
};

class Transcript {
public:
	Transcript() {
		length = 0;
		structure.clear();
		strand = 0;
		seqname = gene_id = transcript_id = "";
		left = "";
	}

	Transcript(const std::string& transcript_id, const std::string& gene_id, const std::string& seqname,
			const char& strand, const std::vector<Interval>& structure, const std::string& left) {
		this->structure = structure;
		this->strand = strand;
		this->seqname = seqname;
		this->gene_id = gene_id;
		this->transcript_id = transcript_id;

		//eliminate prefix spaces in string variable "left"
		int pos = 0;
		int len = left.length();
		while (pos < len && left[pos] == ' ') ++pos;
		this->left = left.substr(pos);

		length = 0;
		int s = structure.size();
		for (int i = 0; i < s; i++) length += structure[i].end + 1 - structure[i].start;
	}

	bool operator< (const Transcript& o) const {
	  return gene_id < o.gene_id || (gene_id == o.gene_id && transcript_id < o.transcript_id) || (gene_id == o.gene_id && transcript_id == o.transcript_id && seqname < o.seqname);
	}

	const std::string& getTranscriptID() const { return transcript_id; }

	const std::string& getGeneID() const { return gene_id; }

	const std::string& getSeqName() const { return seqname; }

	char getStrand() const { return strand; }

	const std::string& getLeft() const { return left; }

	int getLength() const { return length; }

	const std::vector<Interval>& getStructure() const { return structure; }

	void extractSeq (const std::string&, std::string&) const;

	void read(std::ifstream&);
	void write(std::ofstream&);

private:
	int length; // transcript length
	std::vector<Interval> structure; // transcript structure , coordinate starts from 1
	char strand;
	std::string seqname, gene_id, transcript_id; // follow GTF definition
	std::string left;
};

#endif /* TRANSCRIPT_H_ */
