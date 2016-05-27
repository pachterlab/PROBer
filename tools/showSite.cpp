#include<cstdio>
#include<cstring>
#include<cstdlib>
#include<iostream>

#include "genome_structs.h"

using namespace std;

Dictionary dictionary;

int main(int argc, char* argv[]) {
  if (argc != 5) {
    printf("Usage: showSite annotation_name chr strand pos\n");
    exit(-1);
  }

  dictionary.loadIndex(argv[1]);
  cout<< dictionary.showSite(argv[2], argv[3][0], atoi(argv[4]))<< endl;

  return 0;
}
