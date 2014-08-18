#include<ctime>
#include<cstdio>
#include<cstring>
#include<cstdlib>
#include<cassert>
#include<string>
#include<vector>
#include<fstream>

#include "sam/bam.h"
#include "sam/sam.h"

#include "utils.h"
#include "AlignmentGroup.hpp"
#include "SamParser.hpp"
#include "BamWriter.hpp"

using namespace std;

int num_threads;

const bam_header_t *header;
SamParser *parser;
BamWriter *writer;
AlignmentGroup ag;

int main(int argc, char* argv[]) {
  if (argc != 4) { 
    printf("Usage: putTogether input_name num_threads output_name\n");
    exit(-1);
  }

  time_t a = time(NULL);

  num_threads = atoi(argv[2]);

  char inpF[STRLEN], outF[STRLEN];

  sprintf(inpF, "%s_0.bam", argv[1]);
  parser = new SamParser('b', inpF, NULL);
  header = parser->getHeader();
  sprintf(outF, "%s.bam", argv[3]);
  writer = new BamWriter(outF, header, "RSEM");

  while (parser->next(ag)) writer->write(ag, 2);
  printf("Done!\n");

  for (int i = 1; i < num_threads; i++) {
    delete parser;
    sprintf(inpF, "%s_%d.bam", argv[1], i);
    parser = new SamParser('b', inpF, NULL);
    while (parser->next(ag)) writer->write(ag, 2);
    printf("Done %d!\n", i);
  }

  delete parser;
  sprintf(inpF, "%s_N0.bam", argv[1]);
  parser = new SamParser('b', inpF, NULL);
  while (parser->next(ag)) writer->write(ag, 2);
  printf("Done!\n");

  delete parser;
  sprintf(inpF, "%s_N2.bam", argv[1]);
  parser = new SamParser('b', inpF, NULL);
  while (parser->next(ag)) writer->write(ag, 2);
  printf("Done!\n");

  time_t b = time(NULL);
  printf("%d\n", int(b - a));

  return 0;
}
