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
string anno;

ifstream fin;
ofstream fout;

int main(int argc, char* argv[]) {
  if (argc != 4) {
    printf("Usage: showGeneNames annotation_name input.anno.site_info output.gene.site_info\n");
    exit(-1);
  }

  dictionary.loadIndex(argv[1]);

  fin.open(argv[2]);
  fout.open(argv[3]);

  string line;
  istringstream strin;
  int cnt = 0;
  
  while (getline(fin, line)) {
    strin.clear(); strin.str(line);
    strin>> chr>> strand>> pos>> nuniq>> nmulti>> anno;
    fout<< chr<< ' '<< strand<< ' '<< pos<< '\t'<< nuniq<< '\t'<< nmulti<< dictionary.showGeneName(chr, strand, anno)<< endl;

    ++cnt;
    if (cnt % 1000000 == 0) printf("FIN %d\n", cnt);
  }
  fin.close();
  fout.close();
  
  return 0;
}
