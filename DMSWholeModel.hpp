#ifndef DMSWHOLEMODEL_H_
#define DMSWHOLEMODEL_H_

#include<cmath>
#include<cassert>
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
      counts[alignG.aligns[i]->tid] += alignG.aligns[i]->frac; // used for allocating transcripts, frac is 1/total_alignments
    }
  }

  /*
    @comment: Allocate transcripts to threads, calculate auxiliary arrays for each transcript and initialize theta
   */
  void init_for_EM();

  /*
    @param   tid   transcript id
    @param   pos   leftmost position from 5' end, 0-based
    @param   fragment_length  fragment length, 0 means SE read 
    @return   probability of generating such a read (not condition on that read passes the size selection step
   */
  double getProb(int tid, int pos = 0, int fragment_length = 0) const {
    assert(tid >= 0 && tid <= M);
    double prob = theta[tid];
    if (tid > 0) prob *= (fragment_length > 0 ? transcripts[tid]->getProb(pos, fragment_length) : transcripts[tid]->getProb(pos)); 
    return prob;
  }

  /*
    @return   the probability of a read passing the size selection step
   */
  double getProbPass() const {
    return prob_pass;
  }

  /*
    @param   tid   transcript id, 0 means noise transcript
   */
  double getTheta(int tid) {
    assert(tid >= 0 && tid <= M);
    return theta[tid];
  }

  /*
    @param   count0   expected count for backgroud noise
    @param   round    number of EM iterations
   */
  void runEM(double count0, int round);

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
  double prob_pass; // probability of generating a read that pass the size selection step
  std::vector<double> counts, unobserved; // number of observed/unobserved reads fall into each transcript

  double *cdf; // a cumulative array of theta_i * prob_pass_i, used for simulation

  // Params, used for multi-threading
  struct Params {
    int id;
    DMSWholeModel *pointer;

    int num_trans;
    std::vector<DMSTransModel*> trans;
    std::vector<int> origin_ids;

    double *start2, *end2;

    Params(int id, DMSWholeModel *pointer) : id(id), pointer(pointer) {
      num_trans = 0;
      trans.clear();
      origin_ids.clear();
      start2 = end2 = NULL;
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

  /*
    @comment:  Calculate the probability that any read passes the size selection step
   */
  void calcProbPass() {
    prob_pass = theta[0];
    for (int i = 1; i <= M; ++i) 
      if (!isZero(counts[i])) {
	assert(!isZero(theta[i]));
	prob_pass += theta[i] * transcripts[i]->getProbPass();
      }
    assert(!isZero(prob_pass));
  }

  /*
    @comment: This function tries to allocate transcripts to threads evenly
   */
  void allocateTranscriptsToThreads();

  void run_calcAuxiliaryArrays(Params* params) {
    for (int i = 0; i < params->num_trans; ++i)  
      params->trans[i]->calcAuxiliaryArrays();
  }

  void run_makeUpdates(Params* params) {
    for (int i = 0; i < params->num_trans; ++i) 
      params->trans[i]->update();
  }

  void run_EM_step(Params* params) {
    for (int i = 0; i < params->num_trans; ++i) 
      params->trans[i]->EM(N_tot * theta[params->origin_ids[i]]);
  }

  static void* run_calcAuxiliaryArrays_per_thread(void* args) {
    Params *params = (Params*)args;
    params->pointer->run_calcAuxiliaryArrays(params);
    return NULL;
  }

  static void* run_makeUpdates_per_thread(void* args) {
    Params *params = (Params*)args;
    params->pointer->run_makeUpdates(params);
    return NULL;
  }

  static void* run_EM_step_per_thread(void* args) {
    Params *params = (Params*)args;
    params->pointer->run_EM_step(params);
    return NULL;
  }
};

#endif
