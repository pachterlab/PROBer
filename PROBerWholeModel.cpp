#include<cstdio>
#include<cstring>
#include<cassert>
#include<string>
#include<vector>
#include<fstream>
#include<pthread.h>

#include "sam/bam.h"
#include "utils.h"
#include "my_assert.h"
#include "Transcript.hpp"
#include "Transcripts.hpp"
#include "MyHeap.hpp"
#include "PROBerWholeModel.hpp"

PROBerWholeModel::PROBerWholeModel(const char* config_file, int init_state, const Transcripts* trans, int num_threads, int read_length, bool isMAP) {
  // set PROBerTransModel static member values
  int primer_length, min_frag_len, max_frag_len;
  double gamma_init, beta_init;

  FILE *fi = fopen(config_file, "r");
  assert(fi != NULL);
  assert(fscanf(fi, "%d %d %d", &primer_length, &min_frag_len, &max_frag_len) == 3);
  PROBerTransModel::setGlobalParams(primer_length, min_frag_len, max_frag_len, init_state);
  if (trans != NULL) {
    assert(fscanf(fi, "%lf %lf", &gamma_init, &beta_init) == 2);
    PROBerTransModel::setLearningRelatedParams(gamma_init, beta_init, 1.0, read_length, isMAP);
  }
  fclose(fi);

  // initialize data members
  this->num_threads = 0;

  M = 0; 
  theta.clear();
  transcripts.clear();

  N_tot = 0.0;
  threads.clear();

  paramsVecEM.clear();

  for (int i = 0; i < 2; ++i) {
    counts[i].clear();
    unobserved[i].clear();
    prob_noise[i][0] = prob_noise[i][1] = 0.0;
    prob_pass[i] = 0.0;
    paramsVecUp[i].clear();
  }

  sim_tid = -1;
  cdf = NULL;

  channel_to_calc = -1;

  if (trans != NULL) {
    assert(num_threads >= 1);
    this->num_threads = num_threads;

    M = trans->getM();
    theta.assign(M + 1, 0.0);
    transcripts.assign(M + 1, NULL);
    for (int i = 1; i <= M; ++i) {
      const Transcript& tran = trans->getTranscriptAt(i);
      transcripts[i] = new PROBerTransModel(i, tran.getTranscriptID(), tran.getLength());
    }

    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
    threads.assign(num_threads, pthread_t());
    
    int channel = PROBerTransModel::getChannel();
    counts[channel].assign(M + 1, 0.0);
    unobserved[channel].assign(M + 1, 0.0);

    if (PROBerTransModel::isJoint()) {
      counts[channel ^ 1].assign(M + 1, 0.0);
      unobserved[channel ^ 1].assign(M + 1, 0.0);
    }
  }
}

PROBerWholeModel::~PROBerWholeModel() {
  assert(transcripts[0] == NULL);
  for (int i = 1; i <= M; ++i) delete transcripts[i];

  if (PROBerTransModel::isLearning()) pthread_attr_destroy(&attr);

  for (int i = 0; i < 2; ++i)
    for (int j = 0; j < (int)paramsVecUp[i].size(); ++j) delete paramsVecUp[i][j];

  for (int i = 0; i < (int)paramsVecEM.size(); ++i) delete paramsVecEM[i];
}

void PROBerWholeModel::init() {
  int state = PROBerTransModel::getState();
  int channel = PROBerTransModel::getChannel();

  allocateTranscriptsToThreads(state, channel);

  // Initialize theta, prob_noise and auxiliary arrays
  if (state < 3) {
    int total = 1; // noise transcript always counts
    for (int i = 1; i <= M; ++i) 
      if (!transcripts[i]->isExcluded()) ++total;
    prob_noise[channel][0] = 1.0 / total;
    prob_noise[channel][1] = (total - 1.0) / total;
    if (state == 2) 
      memcpy(prob_noise[channel ^ 1], prob_noise[channel], sizeof(double) * 2);
    --total;
    general_assert(total > 0, "No reads aligned to the reference!");
    for (int i = 1; i <= M; ++i) 
      if (!transcripts[i]->isExcluded()) theta[i] = 1.0 / total;

    // run init for each transcript
    int size = paramsVecEM.size();
    channel_to_calc = channel;

    // create threads
    for (int i = 0; i < size; ++i) {
      rc = pthread_create(&threads[i], &attr, run_calcAuxiliaryArrays_per_thread, (void*)paramsVecEM[i]);
      pthread_assert(rc, "pthread_create", "Cannot create thread " + itos(i) + " (numbered from 0) for run_calcAuxiliaryArrays_per_thread at " + cstrtos(channelStr[channel_to_calc]) + " channel!");
    }
    // join threads
    for (int i = 0; i < size; ++i) {
      rc = pthread_join(threads[i], NULL);
      pthread_assert(rc, "pthread_join", "Cannot join thread " + itos(i) + " (numbered from 0) for run_calcAuxiliaryArrays_per_thread at " + cstrtos(channelStr[channel_to_calc]) + " channel!");
    }  

    calcProbPass(channel_to_calc);
  }
}

void PROBerWholeModel::EM_step(double count0) {
  int state = PROBerTransModel::getState();
  int channel = PROBerTransModel::getChannel();
  int size = paramsVecEM.size();
  double N_obs, sum, sum2, value;

  // Update counts
  update(count0);

  // Calculate total number of observed reads
  N_obs = 0.0;
  for (int i = 0; i <= M; ++i) 
    if (!isZero(counts[channel][i])) N_obs += counts[channel][i];
  N_tot = N_obs / prob_pass[channel];

  // Calculate expected hidden reads to each transcript
  for (int i = 0; i <= M; ++i) {
    unobserved[channel][i] = 0.0;
    // If no counts, force the unobserved counts to be 0!
    if (i > 0 && !isZero(counts[channel][i])) unobserved[channel][i] = N_tot * prob_noise[channel][1] * theta[i] * (1.0 - transcripts[i]->getProbPass(channel)); 
  }

  // Estimate new gamma/beta parameters
  // create threads
  for (int i = 0; i < size; ++i) {
    rc = pthread_create(&threads[i], &attr, run_EM_step_per_thread, (void*)(paramsVecEM[i]));
    pthread_assert(rc, "pthread_create", "Cannot create thread " + itos(i) + " (numbered from 0) for run_EM_step_per_thread at " + cstrtos(channelStr[channel]) + " channel!");
  }
  // join threads
  for (int i = 0; i < size; ++i) {
    rc = pthread_join(threads[i], NULL);
    pthread_assert(rc, "pthread_join", "Cannot join thread " + itos(i) + "(numbered from 0) for run_EM_step_per_thread at " + cstrtos(channelStr[channel]) + " channel!");
  }

  // Estimate new theta and prob_noise
  sum = sum2 = 0.0;  
  for (int i = 1; i <= M; ++i) {
    value = counts[channel][i] + unobserved[channel][i];
    sum += value;
    if (state == 3) value += counts[channel ^ 1][i] + unobserved[channel ^ 1][i]; 
    if (state != 2) { theta[i] = value; sum2 += theta[i]; }
  }

  assert(!isZero(sum));
  sum += counts[channel][0];
  prob_noise[channel][0] = counts[channel][0] / sum;
  prob_noise[channel][1] = 1.0 - prob_noise[channel][0];

  if (state != 2) {
    for (int i = 1; i <= M; ++i) theta[i] /= sum2;
  }

  // calculate the probability of a read passing size selection step for next call
  channel_to_calc = (state >= 2 ? (channel ^ 1) : channel); 
  calcProbPass(channel_to_calc);
}

void PROBerWholeModel::wrapItUp(double count0) {
  int state = PROBerTransModel::getState();
  int channel = PROBerTransModel::getChannel();

  // Update counts
  update(count0);

  if (state == 2) {
    int size = paramsVecEM.size();
    channel_to_calc = channel ^ 1;

    // create threads
    for (int i = 0; i < size; ++i) {
      rc = pthread_create(&threads[i], &attr, run_calcAuxiliaryArrays_per_thread, (void*)paramsVecEM[i]);
      pthread_assert(rc, "pthread_create", "Cannot create thread " + itos(i) + " (numbered from 0) for run_calcAuxiliaryArrays_per_thread at " + cstrtos(channelStr[channel_to_calc]) + " channel!");
    }
    // join threads
    for (int i = 0; i < size; ++i) {
      rc = pthread_join(threads[i], NULL);
      pthread_assert(rc, "pthread_join", "Cannot join thread " + itos(i) + " (numbered from 0) for run_calcAuxiliaryArrays_per_thread at " + cstrtos(channelStr[channel_to_calc]) + " channel!");
    }  
    calcProbPass(channel_to_calc);
  }
}

void PROBerWholeModel::read(const char* input_name, const char* statName) {
  char input_param[STRLEN];
  char input_theta[STRLEN];
  std::ifstream fin;

  int state = PROBerTransModel::getState();
  bool learning = PROBerTransModel::isLearning();

  int tmp_M;

  // load gamma
  if (!learning || (learning && state == 1)) {
    sprintf(input_param, "%s.gamma", input_name);
    fin.open(input_param);
    assert(fin.is_open());

    assert(fin>> tmp_M);
    if (!learning) {
      M = tmp_M;
      theta.assign(M + 1, 0.0);
      transcripts.assign(M + 1, NULL);
      for (int i = 1; i <= M; ++i) transcripts[i] = new PROBerTransModel(i);
    }
    else assert(M == tmp_M);

    for (int i = 1; i <= M; ++i) transcripts[i]->read(fin, 0);

    fin.close();
  }

  // load beta
  if (!learning && state == 1) {
    sprintf(input_param, "%s.beta", input_name);
    fin.open(input_param);
    assert(fin.is_open());
    
    assert((fin>> tmp_M) && (tmp_M == M));
    for (int i = 1; i <= M; ++i) transcripts[i]->read(fin, 1);

    fin.close();
  }

  // load theta
  if (!learning) {
    assert(statName != NULL);
    sprintf(input_theta, "%s_%s.theta", statName, channelStr[state & 1]);
    fin.open(input_theta);
    assert(fin.is_open());
      
    assert((fin>> tmp_M) && (tmp_M == M));
    for (int i = 0; i <= M; ++i) assert(fin>> theta[i]);
    fin.close();
  }

  if (verbose) printf("PROBerWholeModel::read is finished!\n");
}

void PROBerWholeModel::writeExprRes(int state, const char* output_name) {
  char exprF[STRLEN];
  double tpm[M + 1], fpkm[M + 1], l_bar;
  
  l_bar = 0.0;
  for (int i = 1; i <= M; ++i) {
    tpm[i] = theta[i] / (transcripts[i]->getLen() + 1.0);
    l_bar += tpm[i];
  }
  assert(l_bar > 0.0);
  l_bar = 1.0 / l_bar;
  for (int i = 1; i <= M; ++i) {
    tpm[i] *= l_bar * 1e6;
    fpkm[i] = tpm[i] * 1e3 / l_bar;
  }

  switch(state) {
  case 0: sprintf(exprF, "%s_minus.expr", output_name); break;
  case 1: sprintf(exprF, "%s_plus.expr", output_name); break;
  case 2: sprintf(exprF, "%s.expr", output_name); break;
  default: assert(false);
  }

  std::ofstream fout(exprF);
  assert(fout.is_open());

  fout.precision(2);
  fout.setf(std::ios::fixed, std::ios::floatfield);

  fout<< "transcript_id\tlength\teffective_length\t"<< (state <= 1 ? "expected_count" : "expected_count_minus\texpected_count_plus")<< "\tTPM\tFPKM"<< std::endl;
  for (int i = 1; i <= M; ++i) {
    fout<< transcripts[i]->getName()<< '\t'<< transcripts[i]->getLen() + transcripts[i]->get_primer_length()<< '\t'<< transcripts[i]->getLen() + 1<< '\t';
    if (state <= 1) fout<< counts[state & 1][i];
    else fout<< counts[0][i]<< '\t'<< counts[1][i];
    fout<< '\t'<< tpm[i]<< '\t'<< fpkm[i]<< '\t'<< std::endl;
  }

  fout.close();
}

void PROBerWholeModel::write(const char* output_name, const char* statName) {
  char output_param[STRLEN];
  char output_theta[STRLEN];
  //  char output_rate[STRLEN];

  std::ofstream fout;

  assert(PROBerTransModel::isLearning());
  int state = PROBerTransModel::getState();
  assert(state < 3); 

  if (state == 0 || state == 2) {
    sprintf(output_param, "%s.gamma", output_name);
    fout.open(output_param);
    assert(fout.is_open());

    fout.precision(10);
    fout.unsetf(std::ios::floatfield);

    fout<< M<< std::endl;
    for (int i = 1; i <= M; ++i) transcripts[i]->write(fout, 0);

    fout.close();

    sprintf(output_theta, "%s_minus.theta", statName);
    fout.open(output_theta);
    assert(fout.is_open());

    fout.precision(10);
    fout.unsetf(std::ios::floatfield);

    fout<< M<< std::endl;
    fout<< prob_noise[0][0];
    for (int i = 1; i <= M; ++i) fout<< '\t'<< prob_noise[0][1] * theta[i];
    fout<< std::endl;
    
    fout.close();
  }

  if (state == 1 || state == 2) {
    sprintf(output_param, "%s.beta", output_name);
    fout.open(output_param);
    assert(fout.is_open());

    fout.precision(10);
    fout.unsetf(std::ios::floatfield);

    fout<< M<< std::endl;
    for (int i = 1; i <= M; ++i) transcripts[i]->write(fout, 1);

    fout.close();

    sprintf(output_theta, "%s_plus.theta", statName);
    fout.open(output_theta);
    assert(fout.is_open());

    fout.precision(10);
    fout.unsetf(std::ios::floatfield);

    fout<< M<< std::endl;
    fout<< prob_noise[1][0];
    for (int i = 1; i <= M; ++i) fout<< '\t'<< prob_noise[1][1] * theta[i];
    fout<< std::endl;
    
    fout.close();
  }

  // write out expression results
  writeExprRes(state, output_name);

  if (verbose) printf("PROBerWholeModel::write is finished!\n");
}

void PROBerWholeModel::startSimulation(int sim_tid) {
  int channel = PROBerTransModel::getChannel();

  this->sim_tid = sim_tid;

  if (sim_tid < 0) {
    cdf = new double[M + 1];
    cdf[0] = theta[0];
    for (int i = 1; i <= M; ++i) {
      cdf[i] = 0.0;
      if (theta[i] > 0.0) {
	assert(transcripts[i]->getEffLen() > 0);
	transcripts[i]->startSimulation();
	cdf[i] = theta[i] * transcripts[i]->getProbPass(channel); 
      }
      cdf[i] += cdf[i - 1];
    }
  }
  else {
    assert(theta[sim_tid] > 0.0 && transcripts[sim_tid]->getEffLen() > 0);
    transcripts[sim_tid]->startSimulation();
  }
}

void PROBerWholeModel::finishSimulation() {
  if (sim_tid < 0) {
    for (int i = 1; i <= M; ++i) 
      if (theta[i] > 0.0) transcripts[i]->finishSimulation();
    delete[] cdf;
  }
  else {
    transcripts[sim_tid]->finishSimulation();
  }
}

void PROBerWholeModel::allocateTranscriptsToThreads(int state, int channel) {
  int id;
  MyHeap my_heap;
  std::vector<int> max_lens; // record maximum len in each thread

  // Allocate transcripts for updating
  my_heap.init(num_threads);
  std::vector<Params*> &paramsVecU = paramsVecUp[channel];
  paramsVecU.assign(num_threads, NULL);
  for (int i = 0; i < num_threads; ++i) paramsVecU[i] = new Params(i, this);

  for (int i = 1; i <= M; ++i) {
    HIT_INT_TYPE numAlign = transcripts[i]->getNumAlignments(channel);
    if (numAlign > 0) {
      id = my_heap.getTop();
      paramsVecU[id]->trans.push_back(transcripts[i]);
      ++paramsVecU[id]->num_trans;
      my_heap.updateTop(numAlign);
    }
  }

  // delete extra paramsVec
  for (id = num_threads - 1; id >= 0 && paramsVecU[id]->num_trans == 0; --id) {
    delete paramsVecU[id];
    paramsVecU[id] = NULL;
  }
  ++id;
  assert(id > 0);
  if (id < num_threads) paramsVecU.resize(id, NULL);

  if (state == 3) return;

  // Allocate transcripts for EM  
  my_heap.init(num_threads);
  paramsVecEM.assign(num_threads, NULL);
  for (int i = 0; i < num_threads; ++i) paramsVecEM[i] = new Params(i, this);
  max_lens.assign(num_threads, 0);

  // allocate transcripts
  for (int i = 1; i <= M; ++i) 
    // If this transcript is not excluded 
    if (!transcripts[i]->isExcluded()) {
      transcripts[i]->init(); // initialize this transcript for learning
      id = my_heap.getTop();
      paramsVecEM[id]->trans.push_back(transcripts[i]);
      ++paramsVecEM[id]->num_trans;
      if (max_lens[id] < transcripts[i]->getLen()) max_lens[id] = transcripts[i]->getLen();
      my_heap.updateTop(transcripts[i]->getLen());
    }

  // delete extra paramsVec
  for (id = num_threads - 1; id >= 0 && paramsVecEM[id]->num_trans == 0; --id) {
    delete paramsVecEM[id];
    paramsVecEM[id] = NULL;
  }
  ++id;
  assert(id > 0);
  if (id < num_threads) paramsVecEM.resize(id, NULL);

  // allocate start2 and end2 for each thread
  for (int i = 0; i < (int)paramsVecEM.size(); ++i) {
    paramsVecEM[i]->start2 = new double[max_lens[i] + 1];
    paramsVecEM[i]->end2 = new double[max_lens[i] + 1];
    for (int j = 0; j < paramsVecEM[i]->num_trans; ++j)
      paramsVecEM[i]->trans[j]->setStart2andEnd2(paramsVecEM[i]->start2, paramsVecEM[i]->end2);
  }
}

void PROBerWholeModel::update(double count0) {
  int channel = PROBerTransModel::getChannel();
  std::vector<Params*> &paramsVecU = paramsVecUp[channel];
  int size = paramsVecU.size();

  // create threads
  for (int i = 0; i < size; ++i) {
    rc = pthread_create(&threads[i], &attr, run_makeUpdates_per_thread, (void*)(paramsVecU[i]));
    pthread_assert(rc, "pthread_create", "Cannot create thread " + itos(i) + " (numbered from 0) for run_makeUpdates_per_thread at " + cstrtos(channelStr[channel]) + " channel!");
  }
  // join threads
  for (int i = 0; i < size; ++i) {
    rc = pthread_join(threads[i], NULL);
    pthread_assert(rc, "pthread_join", "Cannot join thread " + itos(i) + " (numbered from 0) for run_makeUpdates_per_thread at " + cstrtos(channelStr[channel]) + " channel!");
  }    

  // Set counts
  counts[channel][0] = count0;
  for (int i = 1; i <= M; ++i)
    counts[channel][i] = transcripts[i]->getNobs(); // Because N_obs for each channel is stored separately, this operation is safe
}
