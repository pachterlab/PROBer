#include<cstdio>
#include<cstring>
#include<cstdlib>
#include<cassert>
#include<string>
#include<fstream>
#include<sstream>

#include "genome_structs.h"

using namespace std;

Dictionary dictionary;

string chr;
char strand;
int pos, nuniq;
double nmulti, threshold;

ifstream fin;
ofstream fout;

int main(int argc, char* argv[]) {
  if (argc != 5) {
    printf("Usage: annotateSite annotation_name input.site_info threshold output.annotated.site_info\n");
    exit(-1);
  }

  dictionary.loadIndex(argv[1]);
  threshold = atof(argv[3]);

  fin.open(argv[2]);
  fout.open(argv[4]);

  string line;
  istringstream strin;
  int cnt = 0;
  
  while (getline(fin, line)) {
    strin.clear(); strin.str(line);
    strin>> chr>> strand>> pos>> nuniq>> nmulti;
    if (nuniq + nmulti >= threshold) 
      fout<< chr<< ' '<< strand<< ' '<< pos<< '\t'<< nuniq<< '\t'<< nmulti<< '\t'<< dictionary.annotateSite(chr, strand, pos)<< endl;

    ++cnt;
    if (cnt % 1000000 == 0) printf("FIN %d\n", cnt);
  }
  fin.close();
  fout.close();
  
  return 0;
}
