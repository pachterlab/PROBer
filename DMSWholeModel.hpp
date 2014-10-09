#ifndef DMSWHOLEMODEL_H_
#define DMSWHOLEMODEL_H_

#include<cmath>
#include<vector>
#include<pthread.h>

#include "sam/bam.h"
#include "sampling.hpp"
#include "Transcripts.hpp"
#include "InMemoryStructs.hpp"
#include "DMSTransModel.hpp"

class DMSWholeModel {
public:
  /*
    @function   constructor function
    @param   config_file   Configuration file for DMSTransModel static members
    @param   trans         We obtain transcriopt names and lengths from trans
    @param   num_threads   Number of threads allowed to use
   */
  DMSWholeModel(const char* config_file, const Transcripts* trans = NULL, int num_threads = 1);

  /*
    @function   destructor function, release contents of treads and transcripts
   */
  ~DMSWholeModel();

  /*
    @param   alignG   An alignment group, representing a single read's all alignments
   */
  void addAlignments(InMemAlignG& alignG) {
    for (int i = 0; i < alignG.size; ++i) {
      transcripts[alignG.aligns[i]->tid]->addAlignment(alignG.aligns[i]);
      counts[alignG.aligns[i]->tid] += alignG.aligns[i]->frac; // used for allocating transcripts
    }
  }

  /*
    @param   tid       transcript id
    @param   pos       leftmost position in 5' end, 0-based
    @param   frac      the fractional weight of this read
    @comment: This function is for single-end reads
   */
  void update(int tid, int pos, double frac) {
    if (tid > 0) transcripts[tid]->update(pos, frac);
    else counts[tid] += frac;
  }

  /*
    @param   tid   transcript id
    @param   pos   same as the above function
    @param   fragment_length    the estimated fragment length according to the two mates
    @param   frac      same as the above function
    @comment: This function is for paired-end reads
   */
  void update(int tid, int pos, int fragment_length, double frac) {
    if (tid > 0) transcripts[tid]->update(pos, fragment_length, frac);
    else counts[tid] += frac;
  }

  /*
   */
  void init();

  /*
    @comment: This function tries to allocate transcripts to threads evenly
   */
  void allocateTranscriptsToThreads();

  /*
   */
  void runEM(int max_round);

  /*
    @param   input_name   All input files use input_name as their prefixes
   */
  void read(const char* input_name);

  /*
    @param   output_name   All output files use output_name as their prefixes
   */
  void write(const char* output_name);

  /*
    @comment: prepare for simulation
   */
  void startSimulation();

  /*
    @param   sampler           sampler used for simulation
    @param   tid               transcript id
    @param   pos               start position at 5' end
    @param   fragment_length   fragment length 
   */
  void simulate(Sampler* sampler, int& tid, int& pos, int& fragment_length) {
    tid = sampler->sample(cdf, M + 1);
    if (tid == 0) { pos = fragment_length = 0; }
    else transcripts[tid]->simulate(sampler, pos, fragment_length);
  }
    
  /*
    @comment: release memory used for simulation
   */
  void finishSimulation();

private:
  int num_threads; // number of threads we can use
  bool readGamma; // When call read, if we should read gamma
  int M; // Number of transcripts plus noise transcript as tid = 0
  std::vector<double> theta; // M + 1 elements
  std::vector<DMSTransModel*> transcripts; // DMS models for individual transcripts

  double N_tot; // expected total read counts
  std::vector<double> counts; // number of reads fall into each transcript

  // Params, used for multi-threading
  struct Params {
    int id;
    DMSWholeModel *pointer;

    int num_trans;
    std::vector<DMSTransModel*> trans;
    std::vector<int> origin_ids;

    double *start2, *end2;

    double loglik;

    Params(int id, DMSWholeModel *pointer) : id(id), pointer(pointer) {
      num_trans = 0;
      trans.clear();
      origin_ids.clear();
      start2 = end2 = NULL;
      loglik = 0.0;
    }
    
    ~Params() {
      if (start2 != NULL) delete[] start2;
      if (end2 != NULL) delete[] end2;
    }
  };

  std::vector<pthread_t> threads; // pthreads
  pthread_attr_t attr; // pthread attribute
  int rc; // status of pthread running condition
  std::vector<Params*> paramsVec; // parameters used by each thread

  double *cdf; // a cumulative array of theta_i * prob_pass_i, used for simulation

  void run_calcAuxiliaryArrays(Params* params) {
    params->loglik = 0.0;
    for (int i = 0; i < params->num_trans; ++i) { 
      params->trans[i]->calcAuxiliaryArrays();
      params->loglik += params->trans[i]->calcLogLik();
    }
  }

  void run_EM_step(Params* params) {
    params->loglik = 0.0;
    for (int i = 0; i < params->num_trans; ++i) {
      params->trans[i]->EM(N_tot * theta[params->origin_ids[i]]);
      params->loglik += params->trans[i]->calcLogLik();
    }
  }

  static void* run_calcAuxiliaryArrays_per_thread(void* args) {
    Params *params = (Params*)args;
    params->pointer->run_calcAuxiliaryArrays(params);
    return NULL;
  }

  static void* run_EM_step_per_thread(void* args) {
    Params *params = (Params*)args;
    params->pointer->run_EM_step(params);
    return NULL;
  }
};

#endif
