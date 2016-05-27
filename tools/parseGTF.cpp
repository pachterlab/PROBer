#include<cstdio>
#include<cctype>
#include<cstring>
#include<cstdlib>
#include<cassert>
#include<string>
#include<fstream>
#include<map>

#include "GTFItem.h"
#include "genome_structs.h"

using namespace std;

ifstream fin;
ofstream fidx, fout;
GTFItem gtf_line;
string line;

Transcript transcript;
map<BinKeyType, BinValueType> bin_map;
BinValueType *ptr;

// test if we should skip this line
inline bool skip(const string& line) {
  size_t pos = 0, len = line.length();
  while (pos < len && isspace(line[pos])) ++pos;
  return pos >= len || line[pos] == '#'; // skipping if empty line or commented line
}

int main(int argc, char* argv[]) {
  if (argc != 3) {
    printf("Usage: parseGTF input.gtf output_name\n");
    exit(-1);
  }
  
  int cnt = 0;
  
  fin.open(argv[1]);
  ptr = NULL;
  
  while (getline(fin, line)) {
    if (skip(line)) continue;
    gtf_line.parse(line);
    
    if (gtf_line.getFeature() == "exon" || gtf_line.getFeature() == "CDS" || gtf_line.getFeature() == "UTR") {
      gtf_line.parseAttributes(line);

      if (transcript.transcript_id != gtf_line.getTranscriptID()) {
	if (transcript.process()) ptr->add(transcript);
	
	if (transcript.gene_id != gtf_line.getGeneID()) {
	  ptr = &bin_map[BinKeyType(gtf_line.getSeqName(), gtf_line.getStrand())];
	  transcript.resetGene(ptr->get_gid(), gtf_line.getGeneID(), gtf_line.getGeneName(), gtf_line.getSeqName(), gtf_line.getStrand());
	}
	
	transcript.reset(gtf_line.getTranscriptID(), gtf_line.getTranscriptName());
      }

      assert(transcript.chr == gtf_line.getSeqName() && transcript.strand == gtf_line.getStrand());
      transcript.add(gtf_line.getStart(), gtf_line.getEnd(), gtf_line.getFeature());
    }

    ++cnt;
    if (cnt % 100000 == 0) printf("FIN %d\n", cnt);
  }

  if (transcript.process()) ptr->add(transcript);
  
  fin.close();

  printf("Generating output files.\n");
  fidx.open(string(argv[2]) + ".idx");
  fout.open(string(argv[2]) + ".anno");
  
  for (auto&& my_pair : bin_map) {
    my_pair.second.process();
    fidx<< my_pair.first.chr<< '\t'<< my_pair.first.strand<< '\t'<< fout.tellp()<< endl; 
    fout<< my_pair.first.chr<< '\t'<< my_pair.first.strand<< endl;
    my_pair.second.printOut(fout);
  }

  fidx.close();
  fout.close();
  printf("Done.\n");
  
  return 0;
}
