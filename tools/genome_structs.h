#ifndef GENOME_STRUCTS_H_
#define GENOME_STRUCTS_H_

#include<cstdio>
#include<cassert>
#include<string>
#include<fstream>
#include<sstream>
#include<vector>
#include<map>
#include<set>
#include<algorithm>

struct Annotation {
  int gid, tid;
  char type; // 0, exon; 1, intron; 2, CDS; 3, 5'UTR; 4, 3'UTR; 5, UTR.

  Annotation(int gid, int tid, char type) {
    this->gid = gid;
    this->tid = tid;
    this->type = type;
  }
};

struct Unit {
  int start, end; // 0-based [start, end)
  std::vector<Annotation> annotations;

  Unit() {
    start = end = -1;
    annotations.clear();
  }

  Unit(int start, int end, int gid, int tid, char type) {
    this->start = start;
    this->end = end;
    annotations.emplace_back(gid, tid, type);
  }

  Unit(const std::string& line) {
    parse(line);
  }
  
  Unit(const Unit& o) :start(o.start), end(o.end), annotations(o.annotations) {
  }
  
  Unit(Unit&& o) noexcept : start(o.start), end(o.end), annotations(std::move(o.annotations)) {
  }

  Unit& operator= (const Unit& o) {
    start = o.start; end = o.end;
    annotations = o.annotations;
    return *this;
  }

  Unit& operator= (Unit&& o) noexcept {
    start = o.start; end = o.end;
    annotations = std::move(o.annotations);
    return *this;
  }
    
  bool operator< (const Unit& o) const {
    if (start != o.start) return start < o.start;
    return end < o.end;
  }
  
  bool operator== (const Unit& o) const {
    return start == o.start && end == o.end;
  }

  void parse(const std::string& line) {
    int s, gid, tid, type;
    std::istringstream strin(line);

    strin>> start>> end>> s;
    annotations.clear();
    
    for (int i = 0; i < s; ++i) {
      strin>> gid>> tid>> type;
      annotations.emplace_back(gid, tid, (char)type);
    }
  }
  
  std::string toString() {
    std::ostringstream strout;
    strout<< start<< '\t'<< end<< '\t'<< annotations.size();
    for (auto&& annotation : annotations) 
      strout<< '\t'<< annotation.gid<< ' '<< annotation.tid<< ' '<< (int)annotation.type;
    return strout.str();
  }
};

struct Transcript {
  int gid, tid;
  char strand;

  std::string gene_id, gene_name, transcript_id, transcript_name, chr;
  
  std::vector<Unit> exons;
  std::vector<Unit> introns;

  Transcript() {
    gid = tid = -1; strand = 0;
    gene_id = gene_name = transcript_id = transcript_name = chr = "";
    exons.clear(); introns.clear();
  }

  void resetGene(int gid, const std::string& gene_id, const std::string& gene_name, const std::string& chr, char strand) {
    this->gid = gid;
    this->gene_id = gene_id;
    this->gene_name = gene_name;
    this->chr = chr;
    this->strand = strand;

    tid = -1;
  }

  void reset(const std::string& transcript_id, const std::string& transcript_name) {
    ++tid;
    this->transcript_id = transcript_id;
    this->transcript_name = transcript_name;
    exons.clear();
  }

  char get_type(const std::string& feature) {
    if (feature == "exon") return 0;
    if (feature == "CDS") return 2;
    if (feature == "UTR") return 5;
    assert(false);
  }

  void add(int start, int end, const std::string& feature) {
    exons.emplace_back(start - 1, end, gid, tid, get_type(feature));
  }

  // If UTR, exon, CDS or CDS, exon, UTR
  bool is_middle_exon(int i) {
    return exons[i - 1].start == exons[i].start && exons[i].end == exons[i + 1].end; 
  }
  
  bool process() {
    if (tid < 0) return false;

    int s = exons.size();
    int p;
    
    assert(s > 0);
    sort(exons.begin(), exons.end());
    
    // Merge Exons
    p = 0;
    for (int i = 1; i < s; ++i) 
      if (exons[p] == exons[i] || (i + 1 < s && is_middle_exon(i))) {
	assert(exons[p].annotations[0].type == 0 || exons[i].annotations[0].type == 0);
	if (exons[p].annotations[0].type == 0 && exons[i].annotations[0].type > 0)
	  exons[p].annotations[0].type = exons[i].annotations[0].type;
      }
      else {
	++p;
	if (p < i) exons[p] = std::move(exons[i]);
      }
    s = p + 1;
    exons.resize(s);
    
    // modify UTR to be either 3' UTR or 5' UTR
    assert(s > 1 || exons[0].annotations[0].type != 5);
    if (exons[0].annotations[0].type == 5) exons[0].annotations[0].type = (strand == '+' ? 3 : 4);
    if (exons[s - 1].annotations[0].type == 5) exons[s - 1].annotations[0].type = (strand == '+' ? 4 : 3);
    p = 1;
    while (p < s - 1 && exons[p].annotations[0].type == 5)
      exons[p++].annotations[0].type = (strand == '+' ? 3 : 4);
    for (int i = s - 2; i >= p && exons[i].annotations[0].type == 5; --i)
      exons[i].annotations[0].type = (strand == '+' ? 4 : 3);

    // create introns
    introns.clear();
    for (int i = 1; i < s; ++i)
      if (exons[i].start > exons[i - 1].end)
	introns.emplace_back(exons[i - 1].end, exons[i].start, gid, tid, 1);
    
    return true;
  }
};

struct GeneType {
  std::string gene_id, gene_name;
  std::vector<std::string> transcript_ids;
  std::vector<std::string> transcript_names;

  GeneType(const std::string& gene_id, const std::string& gene_name) {
    this->gene_id = gene_id;
    this->gene_name = gene_name;
    transcript_ids.clear();
    transcript_names.clear();
  }

  GeneType(const std::string& line) {
    parse(line);
  }
  
  void parse(const std::string& line) {
    int s;
    std::string value;
    std::istringstream strin(line);
    
    getline(strin, gene_id, '\t');
    getline(strin, gene_name, '\t');
    strin>> s; getline(strin, value, '\t');
    transcript_ids.clear(); transcript_names.clear();
    for (int i = 0; i < s; ++i) {
      getline(strin, value, '\t');
      transcript_ids.push_back(value);
      getline(strin, value, '\t');
      transcript_names.push_back(value);
    }
  }
  
  std::string toString() {
    std::ostringstream strout;
    strout<< gene_id<< '\t'<< gene_name<< '\t'<< transcript_ids.size();
    for (size_t i = 0; i < transcript_ids.size(); ++i)
      strout<< '\t'<< transcript_ids[i]<< '\t'<< transcript_names[i];
    
    return strout.str();
  }
};
  
struct BinKeyType {
  std::string chr;
  char strand;

  BinKeyType() : chr(""), strand(0) {
  }
  
  BinKeyType(const std::string& chr, char strand) : chr(chr), strand(strand) {
  }

  bool operator< (const BinKeyType& o) const {
    if (chr != o.chr) return chr < o.chr;
    return strand < o.strand;
  }
};

struct BinValueType {
  int gid;
  std::vector<GeneType> genes;
  std::vector<Unit> units;
  
  BinValueType() {
    gid = -1;
    genes.clear();
    units.clear();
  }

  int get_gid() { return ++gid; }
  
  void add(Transcript& transcript) {
    if (gid >= genes.size()) genes.emplace_back(transcript.gene_id, transcript.gene_name);
    genes[gid].transcript_ids.push_back(std::move(transcript.transcript_id));
    genes[gid].transcript_names.push_back(std::move(transcript.transcript_name));

    for (auto&& exon : transcript.exons) units.push_back(std::move(exon));
    for (auto&& intron : transcript.introns) units.push_back(std::move(intron));
  }

  void process() {
    int p = 0, s = units.size();
    
    assert(s > 0);
    sort(units.begin(), units.end());

    for (int i = 1; i < s; ++i)
      if (units[p] == units[i])
	for (auto&& anno : units[i].annotations)
	  units[p].annotations.push_back(anno);
      else {
	++p;
	if (p < i) units[p] = std::move(units[i]);
      }

    s = p + 1;
    units.resize(s);
  }

  void loadBin(std::ifstream& fin) {
    int s;
    std::string line;
    
    fin>> s;
    genes.clear(); 
    getline(fin, line);
    for (int i = 0; i < s; ++i) {
      getline(fin, line);
      genes.emplace_back(line);
    }
    fin>> s;
    units.clear();
    getline(fin, line);
    for (int i = 0; i < s; ++i) {
      getline(fin, line);
      units.emplace_back(line);
    }
  }
  
  void printOut(std::ofstream& fout) {
    fout<< gid + 1<< std::endl;
    for (auto&& gene : genes) fout<< gene.toString()<< std::endl;
    fout<< units.size()<< std::endl;
    for (auto&& unit : units) fout<< unit.toString()<< std::endl;
  }
};

struct Entry {
  int uid; // unit id;
  int end; // end position

  Entry(int uid, int end) : uid(uid), end(end) {
  }
  
  bool operator< (const Entry& o) const {
    return end > o.end;
  }
};


struct Dictionary {
  std::map<BinKeyType, long long> idx_map;
  BinKeyType myKey;
  BinValueType myBin;
  std::vector<Entry> my_heap;
  std::ifstream fin_anno;

  bool exist; 
  int p; // a pointer for myBin.units
  
  Dictionary() : exist(false), p(0) {
  }

  ~Dictionary() {
    if (fin_anno.is_open()) fin_anno.close();
  }
  
  void loadIndex(const std::string& file_name) {
    std::string chr;
    char strand;
    long long pos;
    std::ifstream fin(file_name + ".idx");
    
    idx_map.clear();
    while (fin>> chr>> strand>> pos) {
      idx_map[BinKeyType(chr, strand)] = pos;
    }

    fin.close();

    fin_anno.open(file_name + ".anno");
  }

  // 0, no need to load; 1, load fail; 2, load success 
  int loadBin(const std::string& chr, char strand) {
    if (myKey.chr == chr && myKey.strand == strand) return 0;

    myKey.chr = chr; myKey.strand = strand;
    auto it = idx_map.find(myKey);

    if (it == idx_map.end()) return 1;

    std::string tchr;
    char tstrand;
    
    fin_anno.seekg(it->second);
    fin_anno>> tchr>> tstrand;
    assert(tchr == it->first.chr && tstrand == it->first.strand);
    
    myBin.loadBin(fin_anno);
    
    return 2;
  }

  
  std::string annotateSite(const std::string& chr, char strand, int pos) {
    switch(loadBin(chr, strand)) {
    case 0: break;
    case 1: exist = false; break;
    case 2: exist = true; my_heap.clear(); p = 0;
    }

    if (!exist) return "NULL";

    while (!my_heap.empty() && my_heap.front().end <= pos) {
      std::pop_heap(my_heap.begin(), my_heap.end());
      my_heap.pop_back();
    }
    while (p < myBin.units.size() && myBin.units[p].start <= pos) {
      if (myBin.units[p].end > pos) {
	my_heap.emplace_back(p, myBin.units[p].end);
	std::push_heap(my_heap.begin(), my_heap.end());
      }
      ++p;
    }
    
    if (my_heap.empty()) return "None";

    std::string res = "";
    for (auto&& entry : my_heap) {
      if (res != "") res += ",";
      res += std::to_string(entry.uid);
    }

    return res;
  }

  // ShowGeneName provides a string with '\t' in front of each gid
  std::string showGeneName(const std::string& chr, char strand, const std::string& anno) {
    if (anno == "None") return "\tIntergenic";
    if (anno == "NULL") return "\tNo record";

    assert(loadBin(chr, strand) != 1);

    std::string value, res;
    std::istringstream strin(anno);
    std::set<int> gids;
    while (getline(strin, value, ',')) {
      int uid = std::stoi(value);
      for (auto&& anno : myBin.units[uid].annotations) {
	gids.insert(anno.gid);
      }
    }
    res = "";
    for (auto&& gid : gids) {
      res += "\t" + myBin.genes[gid].gene_id;
      if (myBin.genes[gid].gene_name != "") res += "(" + myBin.genes[gid].gene_name + ")";
    }

    return res;
  }
};

#endif
