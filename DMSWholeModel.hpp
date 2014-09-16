#ifndef DMSWHOLEMODEL_H_
#define DMSWHOLEMODEL_H_

#include<vector>

#include "sam/bam.h"
#include "sampling.hpp"
#include "DMSTransModel.hpp"

class DMSWholeModel {
public:
  /*
    @function   constructor function
    @param   config_file   Configuration file for DMSTransModel static members
    @param   header        Bam header, we obtain transcriopt names and lengths from the header
    @param   num_threads   Number of threads allowed to use
   */
  DMSWholeModel(const char* config_file, const bam_header_t* header = NULL, int num_threads = 1);

  /*
    @function   destructor function, release contents of treads and transcripts
   */
  ~DMSWholeModel();

  /*
    @param   tid       transcript id
    @param   pos       leftmost position in 5' end, 0-based
    @param   frac      the fractional weight of this read
    @comment: This function is for single-end reads
   */
  void update(int tid, int pos, double frac) {
    counts[tid] += frac;
    if (tid > 0) transcripts[tid]->update(pos, frac);
  }

  /*
    @param   tid   transcript id
    @param   pos   same as the above function
    @param   fragment_length    the estimated fragment length according to the two mates
    @param   frac      same as the above function
    @comment: This function is for paired-end reads
   */
  void update(int tid, int pos, int fragment_length, double frac) {
    counts[tid] += frac;
    if (tid > 0) transcripts[tid]->update(pos, fragment_length, frac);
  }

  /*
    @comment: set all counters to 0, prepared for update read counts
   */
  void init() {
    counts.assign(M + 1, 0.0);
    for (int i = 1; i <= M; ++i) 
      if (transcripts[i]->canProduceReads()) transcripts[i]->init();
  }

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
    double *start2, *end2;

    double loglik;

    Params(int id, DMSWholeModel *pointer) : id(id), pointer(pointer) {
      num_trans = 0;
      trans.clear();
      start2 = end2 = NULL;
      loglik = 0.0;
    }
    
    ~Params() {
      if (start2 != NULL) delete[] start2;
      if (end2 != NULL) delete[] end2;
    }
  };

  std::vector<Params*> paramsVec; // parameters used by each thread

  double *cdf; // a cumulative array of theta_i * prob_pass_i, used for simulation

  /*
    @comment: This function tries to allocate transcripts to threads evenly
   */
  void allocateTranscriptsToThreads();

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
      params->trans[i]->EM(N_tot * theta[i]);
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
