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
    @param   init_state    the initial state (if (-) or (+) channel, if learn jointly or not)
    @param   trans         We obtain transcript names and lengths from trans
    @param   num_threads   Number of threads allowed to use
    @param   read_length   If set, assuming all reads < read_length are due to adaptor trimming
    @param   isMAP         Use MAP estimates if true
   */
  DMSWholeModel(const char* config_file, int init_state, const Transcripts* trans = NULL, int num_threads = 1, int read_length = -1, bool isMAP = true);

  /*
    @function   destructor function, release contents of treads and transcripts
   */
  ~DMSWholeModel();

  /*
    @comment: change channel
   */
  void flipState() { DMSTransModel::flipState(); }

  /*
    @comment: get channel information from DMSTransModel
   */
  int getChannel() const { return DMSTransModel::getChannel(); }

  /*
    @param   alignG   An alignment group, representing a single read's all alignments
   */
  void addAlignments(InMemAlignG* alignG, InMemAlign* aligns) {
    for (int i = 0; i < alignG->size; ++i) 
      if (!transcripts[aligns[i].tid]->addAlignment(aligns + i)) aligns[i].conprb = -1.0; // This alignment is discarded, mark its conprb as -1.0
  }

  /*
    @param   tid   transcript id
    @param   pos   leftmost position from 5' end, 0-based
    @param   fragment_length  fragment length, 0 means SE read 
    @return   probability of generating such a read (not condition on that read passes the size selection step
   */
  double getProb(int tid, int pos = 0, int fragment_length = 0) const {
    assert(tid >= 0 && tid <= M);
    if (tid == 0) return prob_noise[DMSTransModel::getChannel()][0];
    return prob_noise[DMSTransModel::getChannel()][1] * theta[tid] * (fragment_length > 0 ? transcripts[tid]->getProb(pos, fragment_length) : transcripts[tid]->getProb(pos)); 
  }

  /*
    @return   the probability of a read passing the size selection step
   */
  double getProbPass() const {
    return prob_pass[DMSTransModel::getChannel()];
  }

  /*
    @param   tid   transcript id, 0 means noise transcript
    @return   the theta defined as the fraction of reads, including noise transcript
   */
  double getTheta(int tid) {
    assert(tid >= 0 && tid <= M);
    if (tid == 0) return prob_noise[DMSTransModel::getChannel()][0];
    return prob_noise[DMSTransModel::getChannel()][1] * theta[tid];
  }

  /*
    @comment: Allocate transcripts to threads, calculate auxiliary arrays by calling DMSTransModel::init() for each transcript and initialize theta
   */
  void init();

  /*
    @param   count0   expected count for backgroud noise
    @comment: Run one iteration of EM algorithm on the transcriptome 
   */
  void EM_step(double count0);

  /*
    @param   count0   expected count for backgroud noise
    @comment: Run one iteration of EM algorithm on the transcriptome 
   */
  void wrapItUp(double count0);

  /*
    @param   input_name   All input files use input_name as their prefixes
   */
  void read(const char* input_name);

  /*
    @param   output_name   All output files use output_name as their prefixes
   */
  void write(const char* output_name);

  /*
    @param   sim_tid   if only simulate reads from sim_tid, default is not (-1)
    @comment: prepare for simulation
   */
  void startSimulation(int sim_tid = -1);

  /*
    @param   sampler           sampler used for simulation
    @param   tid               transcript id
    @param   pos               start position at 5' end
    @param   fragment_length   fragment length 
   */
  void simulate(Sampler* sampler, int& tid, int& pos, int& fragment_length) {
    if (sim_tid < 0) {
      tid = sampler->sample(cdf, M + 1);
      if (tid == 0) { pos = fragment_length = 0; }
      else transcripts[tid]->simulate(sampler, pos, fragment_length);
    }
    else {
      tid = sim_tid;
      transcripts[tid]->simulate(sampler, pos, fragment_length);
    }
  }
    
  /*
    @comment: release memory used for simulation
   */
  void finishSimulation();

private:
  int num_threads; // number of threads we can use
  int M; // Number of transcripts
  std::vector<double> theta; // M + 1 elements. If we learn parameters, theta[0] = 0 and sum_{i=1}^{M} theta[i] = 1, the actual fraction of i is prob_noise[channel][i] * theta[i]; However, if we simulate, theta[0] > 0, theta[i] represents the actual fraction of reads coming from transcript i and sum_{i=0}^{M} theta[i] = 1.
  std::vector<DMSTransModel*> transcripts; // DMS models for individual transcripts

  double N_tot; // expected total read counts
  std::vector<double> counts[2], unobserved[2]; // number of observed/unobserved reads fall into each transcript for two channels

  double prob_noise[2][2]; // the first dimension represent state, the second dimension: 0, probability of generating a noise read; 1, probability of generating a read from transcripts.
  double prob_pass[2]; // probability of generating a read that pass the size selection step

  int sim_tid; // if sim_tid > 0, only simulate from transcript sim_tid
  double *cdf; // a cumulative array of theta_i * prob_pass_i, used for simulation


  int channel_to_calc; // the channel to calculate auxiliary arrays


  // Params, used for multi-threading
  struct Params {
    int id;
    DMSWholeModel *pointer;

    int num_trans;
    std::vector<DMSTransModel*> trans;

    double *start2, *end2;

    Params(int id, DMSWholeModel *pointer) : id(id), pointer(pointer) {
      num_trans = 0;
      trans.clear();
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
  std::vector<Params*> paramsVecUp[2], paramsVecEM; // parameters used by each thread for updates and EM steps

  /*
    @param   channel   (+) or (-), which channel we are working on
    @comment:  Calculate the probability that any read passes the size selection step
   */
  void calcProbPass(int channel) {
    prob_pass[channel] = prob_noise[channel][0];
    for (int i = 1; i <= M; ++i) 
      prob_pass[channel] += prob_noise[channel][1] * theta[i] * transcripts[i]->getProbPass(channel);
    assert(!isZero(prob_pass[channel]));
  }

  /*
    @param   state     state to deal with
    @param   channel   channel to deal with
    @comment: This function tries to allocate transcripts to threads evenly
   */
  void allocateTranscriptsToThreads(int state, int channel);

  /*
    @param   count0    expected count for backgroud noise
   */
  void update(double count0);

  /*
    @param  state   the current state
    @param  output_name   the output name prefix for this data set
    @comment: this procedure write out a table with column names as transcript_id, length, effective_length, expected_count (or expected_count_plus and expected_count_minus), TPM and FPKM
   */
  void writeExprRes(int state, const char* output_name);

  void run_calcAuxiliaryArrays(Params* params) {
    for (int i = 0; i < params->num_trans; ++i)  
      params->trans[i]->calcAuxiliaryArrays(channel_to_calc);
  }

  void run_makeUpdates(Params* params) {
    for (int i = 0; i < params->num_trans; ++i) 
      params->trans[i]->update();
  }

  void run_EM_step(Params* params) {
    int channel = DMSTransModel::getChannel();
    for (int i = 0; i < params->num_trans; ++i) 
      params->trans[i]->EM_step(N_tot * prob_noise[channel][1] * theta[params->trans[i]->getTid()]);
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
