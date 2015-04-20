#include<cmath>
#include<cstdio>
#include<cstring>
#include<cstdlib>
#include<cassert>
#include<string>
#include<vector>
#include<set>
#include<fstream>
#include<pthread.h>

#include "sam/sam.h"
#include "sam/bam.h"

#include "utils.h"
#include "my_assert.h"
#include "MyHeap.hpp"
#include "InMemoryStructs.hpp"
#include "PROBerTransModelS.hpp"

using namespace std;

struct ATranscript {
  PROBerTransModelS *model;
  double Nobs[2]; // minus and plus channel observed reads

  ATranscript() {
    model = NULL;
    Nobs[0] = Nobs[1] = 0.0;
  }

  ~ATranscript() {
    if (model != NULL) delete model;
  }
};

struct TransVecPerThread {
  int no; // thread id
  int s; // total number of transcripts in this thread
  vector<ATranscript*> trans; // transcript vector

  TransVecPerThread() {
    no = -1; s = 0; trans.clear();
  }
};


bool isJoint;

int read_length;
int rounds;
bool paired_end;
int nthreads;
char inputList[STRLEN];
bool isMAP;

int M; // total number of transcripts
vector<ATranscript*> trans; // link to transcript model of each transcript, can be null
 
// thread-related
int rc;
pthread_attr_t attr;
vector<pthread_t> threads;
vector<TransVecPerThread> transvec;

void setupConfig(char* configF) {
  FILE *fi;
  int primer_length, min_frag_len, max_frag_len;
  double gamma_init, beta_init;

  fi = fopen(configF, "r");
  assert(fi != NULL);
  assert(fscanf(fi, "%d %d %d %lf %lf", &primer_length, &min_frag_len, &max_frag_len, &gamma_init, &beta_init) == 5);
  fclose(fi);

  assert(isJoint && isMAP);
  PROBerTransModelS::setGlobalParams(primer_length, min_frag_len, max_frag_len, 2);
  PROBerTransModelS::setLearningRelatedParams(gamma_init, beta_init, 1.0, read_length, isMAP);
}

void assignTranscripts(char* inpBamF, char* listF) {
  samfile_t *in;
  bam_header_t *header;
  set<string> inList;
  MyHeap aheap;
  ATranscript *atran;

  // Load list of transcripts that we are interested, if listF == NULL, all transcripts are considered
  inList.clear();
  if (listF[0] != 0) {
    ifstream fin(listF);
    string key;
    while (getline(fin, key)) {
      inList.insert(key);
    }
    fin.close();
  }

  in = samopen(inpBamF, "rb", NULL);
  header = in->header;

  // Initialize parallel computing required data structures
  threads.assign(nthreads, pthread_t());
  transvec.assign(nthreads, TransVecPerThread());
  for (int i = 0; i < nthreads; ++i) transvec[i].no = i;

  M = header->n_targets;
  trans.assign(M, NULL);
  aheap.init(nthreads);

  int pid;
  for (int i = 0; i < M; ++i)
    if (listF[0] == 0 || inList.find(string(header->target_name[i])) != inList.end()) {
      atran = new ATranscript();
      atran->model = new PROBerTransModelS(i, header->target_name[i], header->target_len[i]);
      trans[i] = atran;
      pid = aheap.getTop();
      aheap.updateTop(atran->model->getLen());
      transvec[pid].trans.push_back(atran);
      ++transvec[pid].s;
    }

  samclose(in);
}

void parseAlignments(char *inpF, int channel) {
  samfile_t *in;
  bam1_t *b;
  InMemAlign ima;
  bool is_paired;

  int cnt = 0;

  in = samopen(inpF, "rb", NULL);
  b = bam_init1();

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

    if (trans[b->core.tid] == NULL) continue;
    
    ima.pos = b->core.pos;
    ima.fragment_length = (is_paired ? abs(b->core.isize) : (b->core.l_qseq < read_length ? b->core.l_qseq : 0));
    assert(ima.fragment_length >= 0);
    uint8_t *p_tag = bam_aux_get(b, "ZW");
    assert(p_tag != NULL);
    ima.frac = double(bam_aux2f(p_tag));
    
    trans[b->core.tid]->model->addAlignment(&ima);
    trans[b->core.tid]->Nobs[channel] += ima.frac;
  }

  // flip state for the next step
  for (int i = 0; i < M; ++i) 
    if (trans[i] != NULL) trans[i]->model->flipState();

  bam_destroy1(b);
  samclose(in);
}

void* runEM(void* arg) {
  TransVecPerThread *params = (TransVecPerThread*)arg;
  ATranscript *atran;
  PROBerTransModelS *model;
  int channel;

  for (int i = 0; i < params->s; ++i) {
    atran = params->trans[i];
    model = atran->model;

    // initialize
    model->init();
  
    // joint mode, EM
    channel = model->getChannel();
    model->calcAuxiliaryArrays(channel);
    for (int ROUND = 0; ROUND < rounds; ++ROUND) {
      model->EM_step(atran->Nobs[channel] / model->getProbPass(channel));
      model->flipState();
      channel = model->getChannel();
      model->EM_step(atran->Nobs[channel] / model->getProbPass(channel));
      model->flipState();
      channel = model->getChannel();
    }

    if ((i + 1) % 50 == 0) printf("%d transcripts are processed in thread %d!\n", (i + 1), params->no);
  }
  
  return NULL;
}

void writeItOut(char* outName) {
  char outF[STRLEN];
  ofstream fout;

  sprintf(outF, "%s.gamma", outName);
  fout.open(outF);
  fout.precision(10);
  fout.unsetf(std::ios::floatfield);
  fout<< M<< endl;
  for (int i = 0; i < M; ++i) {
    if (trans[i] != NULL) trans[i]->model->write(fout, 0);
    else fout<< endl;
  }
  fout.close();
  printf("Gamma file is written!\n");

  sprintf(outF, "%s.beta", outName);
  fout.open(outF);
  fout.precision(10);
  fout.unsetf(std::ios::floatfield);
  fout<< M<< endl;
  for (int i = 0; i < M; ++i) {
    if (trans[i] != NULL) trans[i]->model->write(fout, 1);
    else fout<< endl;
  }
  fout.close();
  printf("Beta file is written!\n");
}

void release() {
  for (int i = 0; i < M; ++i) {
    if (trans[i] != NULL) delete trans[i];
  }
}

int main(int argc, char* argv[]) {
  if (argc < 5) {
    printf("Usage: PROBer_single_transcript_batch config_file minus_channel.bam plus_channel.bam output_name [--read-length read_length] [--rounds number_of_rounds] [--paired-end] [--input input_list.txt] [-p number_of_threads] [--maximum-likelihood]\n");
    exit(-1);
  }

  isJoint = true;

  read_length = -1;
  rounds = 400;
  paired_end = false;
  inputList[0] = 0;
  nthreads = 1;
  isMAP = true;

  for (int i = 5; i < argc; ++i) {
    if (!strcmp(argv[i], "--read-length")) {
      read_length = atoi(argv[i + 1]);
    }
    if (!strcmp(argv[i], "--rounds")) {
      rounds = atoi(argv[i + 1]);
    }
    if (!strcmp(argv[i], "--paired-end")) {
      paired_end = true;
    }
    if (!strcmp(argv[i], "--input")) {
      strcpy(inputList, argv[i + 1]);
    }
    if (!strcmp(argv[i], "-p")) {
      nthreads = atoi(argv[i + 1]);
    }
    if (!strcmp(argv[i], "--maximum-likelihood")) {
      isMAP = false;
    }
  }

  // set up global parameters
  setupConfig(argv[1]);
  printf("Set up config is done!\n");

  // assign transcripts to different threads
  assignTranscripts(argv[2], inputList);
  printf("Assign transcripts is done!\n");

  // parse alignments
  parseAlignments(argv[2], 0);
  printf("parseAlignments for '-' is done!\n");
  parseAlignments(argv[3], 1);
  printf("parseAlignments for '+' is done!\n");

  // Run EM
  for (int i = 0; i < nthreads; ++i) {
    rc = pthread_create(&threads[i], &attr, runEM, (void*)(&transvec[i]));
    pthread_assert(rc, "pthread_create", "Cannot create thread " + itos(i) + " (numbered from 0)!");
  }
  for (int i = 0; i < nthreads; ++i) {
    rc = pthread_join(threads[i], NULL);
    pthread_assert(rc, "pthread_join", "Cannot join thread " + itos(i) + " (numbered from 0)!");
  }
  printf("Run EM is done!\n");

  // output
  writeItOut(argv[4]);
  printf("Results are written!\n");

  // release resource
  release();

  return 0;
}
