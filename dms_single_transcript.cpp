#include<cmath>
#include<cstdio>
#include<cstring>
#include<cstdlib>
#include<cassert>
#include<string>
#include<vector>
#include<fstream>

#include "sam/sam.h"
#include "sam/bam.h"

#include "utils.h"
#include "InMemoryStructs.hpp"
#include "DMSTransModelS.hpp"

using namespace std;

DMSTransModelS *model;

char inF[STRLEN], outF[STRLEN];

FILE *fi, *fo;
ifstream fin;
ofstream fout;

int primer_length, min_frag_len, max_frag_len;
double gamma_init, beta_init;

int read_length;
bool isMAP;
int rounds;
bool isJoint;

int state, channel;

double N_obs[2];

bool paired_end;

inline void addAlignments(DMSTransModelS *model, int channel, string& readname, vector<InMemAlign*>& imas) {
  int size;
  size = imas.size();
  if (size > 1) {
    fprintf(stderr, "Read %s aligns to %d positions in the transcript!\n", readname.c_str(), size);
    for (int i = 0; i < size; ++i) imas[i]->frac = 1.0 / size;
  }
  for (int i = 0; i < size; ++i) {   
    model->addAlignment(imas[i]);
    N_obs[channel] += imas[i]->frac;
  }
}

void loadAlignments(char *inpF, int channel) {
  samfile_t *in;
  bam_header_t *header;
  bam1_t *b;
  InMemAlign *ima;
  vector<InMemAlign*> imas;
  string cname, qname;
  bool is_paired;

  in = samopen(inpF, "rb", NULL);
  header = in->header;
  assert(header->n_targets == 1);
  if (model == NULL) {
    model = new DMSTransModelS(0, header->target_name[0], header->target_len[0]);
  }
  
  b = bam_init1();
  imas.clear();
  cname = "";

  int cnt = 0;

  while (samread(in, b) >= 0) {
    ++cnt;
    if (cnt % 1000000 == 0) printf("Loaded %d lines for %s channel!\n", cnt, channelStr[channel]);
    
    if (b->core.flag & 0x0004) continue; // If unmapped, continue
    is_paired = (b->core.flag & 0x0001);
    
    if (is_paired != paired_end) {
      fprintf(stderr, "Detected a %s read in %s data!\n", (is_paired ? "paired-end" : "single-end"), (paired_end ? "paired-end" : "single-end"));
      exit(-1);
    }

    if (is_paired && !(b->core.flag & 0x0040)) continue; // If paired-end and not the first mate, continue
    assert(!(b->core.flag & 0x0010));
    //    if (b->core.flag & 0x0010) continue; // If read aligns to the reverse strand

    qname = string(bam1_qname(b));
    if (cname != qname) {
      addAlignments(model, channel, cname, imas);
      cname = qname;
      imas.clear();
    }
    ima = new InMemAlign();
    ima->pos = b->core.pos;
    ima->fragment_length = (is_paired ? abs(b->core.isize) : (b->core.l_qseq < read_length ? b->core.l_qseq : 0));
    assert(ima->fragment_length >= 0);
    ima->frac = 1.0;
    imas.push_back(ima);
  }
  addAlignments(model, channel, cname, imas);

  bam_destroy1(b);
  samclose(in);
}

int main(int argc, char* argv[]) {
  if (argc < 5) {
    printf("Usage: dms_single_transcript config_file minus_channel.bam plus_channel.bam output_name [--MAP] [--rounds number_of_rounds] [--read-length read_length] [--sep] [--paired-end]\n");
    exit(-1);
  }

  read_length = -1;
  isMAP = false;
  rounds = 400;
  isJoint = true;
  paired_end = false;
  for (int i = 5; i < argc; ++i) {
    if (!strcmp(argv[i], "--MAP")) isMAP = true;
    if (!strcmp(argv[i], "--rounds")) {
      rounds = atoi(argv[i + 1]);
    }
    if (!strcmp(argv[i], "--read-length")) {
      read_length = atoi(argv[i + 1]);
    }
    if (!strcmp(argv[i], "--sep")) {
      isJoint = false;
    }
    if (!strcmp(argv[i], "--paired-end")) {
      paired_end = true;
    }
  }

  fi = fopen(argv[1], "r");
  assert(fi != NULL);
  assert(fscanf(fi, "%d %d %d %lf %lf", &primer_length, &min_frag_len, &max_frag_len, &gamma_init, &beta_init) == 5);
  fclose(fi);

  DMSTransModelS::setGlobalParams(primer_length, min_frag_len, max_frag_len, (isJoint ? 2 : 0));
  DMSTransModelS::setLearningRelatedParams(gamma_init, beta_init, 1.0, read_length, isMAP);

  model = NULL;
  N_obs[0] = N_obs[1] = 0.0;

  if (isJoint) {
    // load alignments
    loadAlignments(argv[2], DMSTransModelS::getChannel());
    DMSTransModelS::flipState();
    loadAlignments(argv[3], DMSTransModelS::getChannel());
    DMSTransModelS::flipState();

    // initialize
    model->init();
    
    // joint mode, EM
    channel = DMSTransModelS::getChannel();
    model->calcAuxiliaryArrays(channel);
    for (int i = 0; i < rounds; ++i) {
      model->EM_step(N_obs[channel] / model->getProbPass(channel));
      DMSTransModelS::flipState();
      channel = DMSTransModelS::getChannel();
      model->EM_step(N_obs[channel] / model->getProbPass(channel));
      DMSTransModelS::flipState();
      channel = DMSTransModelS::getChannel();
      if ((i + 1) % 100 == 0) printf("FIN %d iterations.\n", i + 1);
    }

    // output
    sprintf(outF, "%s.gamma", argv[4]);
    fout.open(outF);
    fout.precision(10);
    fout.unsetf(std::ios::floatfield);
    model->write(fout, 0);
    fout.close();

    sprintf(outF, "%s.beta", argv[4]);
    fout.open(outF);
    fout.precision(10);
    fout.unsetf(std::ios::floatfield);
    model->write(fout, 1);
    fout.close();

    // release resource
    delete model;
  }
  else {
    // load alignments from minus channel
    channel = DMSTransModelS::getChannel();
    assert(channel == 0);
    loadAlignments(argv[2], channel);

    // initialize
    model->init();

    // EM for (-)
    model->calcAuxiliaryArrays(channel);
    for (int i = 0; i < rounds; ++i) {
      model->EM_step(N_obs[channel] / model->getProbPass(channel));
      if ((i + 1) % 100 == 0) printf("FIN %d iterations for %s channel.\n", i + 1, channelStr[channel]);
    }

    // output gamma
    sprintf(outF, "%s.gamma", argv[4]);
    fout.open(outF);
    fout.precision(10);
    fout.unsetf(std::ios::floatfield);
    model->write(fout, channel);
    fout.close();

    // release resource
    delete model;
    model = NULL;

    // switch to plus channel
    DMSTransModelS::flipState();
    channel = DMSTransModelS::getChannel();
    assert(channel == 1);

    // load alignments from plus channel and new model
    loadAlignments(argv[3], channel);
   
    // load estimated gamma values
    sprintf(inF, "%s.gamma", argv[4]);
    fin.open(inF);
    model->read(fin, 0); // load gammas
    fin.close();

    // initialize
    model->init();

    // EM for (+)
    model->calcAuxiliaryArrays(channel);
    for (int i = 0; i < rounds; ++i) {
      model->EM_step(N_obs[channel] / model->getProbPass(channel));
      if ((i + 1) % 100 == 0) printf("FIN %d iterations for %s channel\n", i + 1, channelStr[channel]);
    }

    // output beta
    sprintf(outF, "%s.beta", argv[4]);
    fout.open(outF);
    fout.precision(10);
    fout.unsetf(std::ios::floatfield);
    model->write(fout, channel);
    fout.close();

    // release resource
    delete model;
  }

  return 0;
}
