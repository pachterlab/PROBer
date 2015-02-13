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
#include "DMSWholeModel.hpp"

DMSWholeModel::DMSWholeModel(const char* config_file, int init_state, const Transcripts* trans, int num_threads, int read_length, bool isMAP) {
  // set DMSTransModel static member values
  int primer_length, min_frag_len, max_frag_len;
  double gamma_init, beta_init;

  FILE *fi = fopen(config_file, "r");
  assert(fi != NULL);
  assert(fscanf(fi, "%d %d %d", &primer_length, &min_frag_len, &max_frag_len) == 3);
  DMSTransModel::setGlobalParams(primer_length, min_frag_len, max_frag_len, init_state);
  if (trans != NULL) {
    assert(fscanf(fi, "%lf %lf", &gamma_init, &beta_init) == 2);
    DMSTransModel::setLearningRelatedParams(gamma_init, beta_init, 1.0, read_length, isMAP);
  }
  fclose(fi);

  // initialize data members
  this->num_threads = 0;

  M = 0; 
  theta.clear();
  transcripts.clear();

  N_tot = 0.0;
  threads.clear();

  for (int i = 0; i < 2; ++i) {
    counts[i].clear();
    unobserved[i].clear();
    prob_noise[i][0] = prob_noise[i][1] = 0.0;
    prob_pass[i] = 0.0;
    paramsVecUp[i].clear();
    paramsVecEM[i].clear();
  }

  cdf = NULL;
  
  if (trans != NULL) {
    assert(num_threads >= 1);
    this->num_threads = num_threads;

    M = trans->getM();
    theta.assign(M + 1, 0.0);
    transcripts.assign(M + 1, NULL);
    for (int i = 1; i <= M; ++i) {
      const Transcript& tran = trans->getTranscriptAt(i);
      transcripts[i] = new DMSTransModel(i, tran.getTranscriptID(), tran.getLength());
    }

    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
    threads.assign(num_threads, pthread_t());
    
    for (int i = 1; i <= M; ++i) theta[i] = 1.0 / M;

    int channel = DMSTransModel::getChannel();
    counts[channel].assign(M + 1, 0.0);
    unobserved[channel].assign(M + 1, 0.0);
    prob_noise[channel][0] = 1.0 / (M + 1);
    prob_noise[channel][1] = M * 1.0 / (M + 1);

    if (DMSTransModel::isJoint()) {
      counts[channel ^ 1].assign(M + 1, 0.0);
      unobserved[channel ^ 1].assign(M + 1, 0.0);
      prob_noise[channel ^ 1][0] = 1.0 / (M + 1);
      prob_noise[channel ^ 1][1] = M * 1.0 / (M + 1);
    }
  }
}

DMSWholeModel::~DMSWholeModel() {
  assert(transcripts[0] == NULL);
  for (int i = 1; i <= M; ++i) delete transcripts[i];

  if (DMSTransModel::isLearning()) pthread_attr_destroy(&attr);

  for (int i = 0; i < 2; ++i) {
    for (int j = 0; j < (int)paramsVecUp[i].size(); ++j) delete paramsVecUp[i][j];
    for (int j = 0; j < (int)paramsVecEM[i].size(); ++j) delete paramsVecEM[i][j];
  }
}

void DMSWholeModel::init() {
  int channel = DMSTransModel::getChannel();
  allocateTranscriptsToThreads(channel);

  // run init for each transcript
  std::vector<Params*> &paramsVecE = paramsVecEM[channel];
  // create threads
  for (int i = 0; i < num_threads; ++i) {
    rc = pthread_create(&threads[i], &attr, run_init_per_thread, (void*)paramsVecE[i]);
    pthread_assert(rc, "pthread_create", "Cannot create thread " + itos(i) + " (numbered from 0) for run_init_per_thread!");
  }
  // join threads
  for (int i = 0; i < num_threads; ++i) {
    rc = pthread_join(threads[i], NULL);
    pthread_assert(rc, "pthread_join", "Cannot join thread " + itos(i) + " (numbered from 0) for run_init_per_thread!");
  }  

  calcProbPass(channel);
}

void DMSWholeModel::update(double count0) {
  int channel = DMSTransModel::getChannel();
  std::vector<Params*> &paramsVecU = paramsVecUp[channel];

  // create threads
  for (int i = 0; i < (int)paramsVecU.size(); ++i) {
    rc = pthread_create(&threads[i], &attr, run_makeUpdates_per_thread, (void*)(paramsVecU[i]));
    pthread_assert(rc, "pthread_create", "Cannot create thread " + itos(i) + " (numbered from 0) for run_makeUpdates_per_thread!");
  }
  // join threads
  for (int i = 0; i < (int)paramsVecU.size(); ++i) {
    rc = pthread_join(threads[i], NULL);
    pthread_assert(rc, "pthread_join", "Cannot join thread " + itos(i) + " (numbered from 0) for run_makeUpdates_per_thread!");
  }    

  // Set counts
  counts[channel][0] = count0;
  for (int i = 1; i <= M; ++i)
    counts[channel][i] = transcripts[i]->getNobs(); // Because N_obs for each channel is stored separately, this operation is safe
}

void DMSWholeModel::EM_step(double count0) {
  int state = DMSTransModel::getState();
  int channel = DMSTransModel::getChannel();
  std::vector<Params*> &paramsVecE = paramsVecEM[channel];
  double N_obs, sum, sum2, value;

  // Update counts
  update(count0);

  // Calculate total number of observed reads
  N_obs = 0.0;
  for (int i = 0; i <= M; ++i) 
    if (isZero(counts[channel][i])) counts[channel][i] = 0.0;
    else {
      N_obs += counts[channel][i];
    }
  N_tot = N_obs / prob_pass[channel];

  // Calculate expected hidden reads to each transcript
  for (int i = 0; i <= M; ++i) {
    unobserved[channel][i] = 0.0;
    // If no counts, force the unobserved counts to be 0!
    if (i > 0 && !isZero(counts[channel][i])) unobserved[channel][i] = N_tot * prob_noise[channel][1] * theta[i] * (1.0 - transcripts[i]->getProbPass()); 
  }
  
  // Estimate new gamma/beta parameters
  // create threads
  for (int i = 0; i < num_threads; ++i) {
    rc = pthread_create(&threads[i], &attr, run_EM_step_per_thread, (void*)(paramsVecE[i]));
    pthread_assert(rc, "pthread_create", "Cannot create thread " + itos(i) + " (numbered from 0) for run_EM_step_per_thread!");
  }
  // join threads
  for (int i = 0; i < num_threads; ++i) {
    rc = pthread_join(threads[i], NULL);
    pthread_assert(rc, "pthread_join", "Cannot join thread " + itos(i) + "(numbered from 0) for run_EM_step_per_thread!");
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

  if (state != 2)
    for (int i = 1; i <= M; ++i) theta[i] /= sum2;

  // calculate the probability of a read passing size selection step for next call
  calcProbPass(DMSTransModel::isJoint() ? channel ^ 1 : channel);
}

void DMSWholeModel::read(const char* input_name) {
  char input_param[STRLEN];
  char input_theta[STRLEN];
  std::ifstream fin;

  int state = DMSTransModel::getState();
  bool learning = DMSTransModel::isLearning();

  int tmp_M;

  if ((state == 0 && !learning) || (state == 1 && learning)) {
    sprintf(input_param, "%s.gamma", input_name);
    fin.open(input_param);
    assert(fin.is_open());

    assert(fin>> tmp_M);
    if (!learning) {
      M = tmp_M;
      theta.assign(M + 1, 0.0);
      transcripts.assign(M + 1, NULL);
      for (int i = 1; i <= M; ++i) transcripts[i] = new DMSTransModel(i);
    }
    else assert(M == tmp_M);

    for (int i = 1; i <= M; ++i) transcripts[i]->read(fin);

    fin.close();

    if (!learning) {
      // Loading expression of minus chanel for the sake of simulation 
      sprintf(input_theta, "%s_minus.theta", input_name);
      fin.open(input_theta);
      assert(fin.is_open());
      
      assert((fin>> tmp_M) && (tmp_M == M));
      for (int i = 0; i <= M; ++i) assert(fin>> theta[i]);
      fin.close();
    }
  }
  else {
    sprintf(input_param, "%s.beta", input_name);
    fin.open(input_param);
    assert(fin.is_open());
    
    assert((fin>> tmp_M) && (tmp_M == M));
    for (int i = 1; i <= M; ++i) transcripts[i]->read(fin);

    fin.close();

    sprintf(input_theta, "%s_plus.theta", input_name);
    fin.open(input_theta);
    assert(fin.is_open());

    assert((fin>> tmp_M) && (tmp_M == M));
    for (int i = 0; i <= M; ++i) assert(fin>> theta[i]);
    fin.close();
  }

  if (verbose) printf("DMSWHoleModel::read is finished!\n");
}

void DMSWholeModel::writeExprRes(int state, const char* output_name) {
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
  case 2: assert(false);
  case 3: sprintf(exprF, "%s.expr", output_name); break;
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

void DMSWholeModel::write(const char* output_name) {
  char output_param[STRLEN];
  char output_theta[STRLEN];
  char output_rate[STRLEN];

  std::ofstream fout;

  assert(DMSTransModel::isLearning());
  int state = DMSTransModel::getState();
  assert(state != 2); 

  if (state == 0 || state == 3) {
    sprintf(output_param, "%s.gamma", output_name);
    fout.open(output_param);
    assert(fout.is_open());

    fout.precision(10);
    fout.unsetf(std::ios::floatfield);

    fout<< M<< std::endl;
    for (int i = 1; i <= M; ++i) transcripts[i]->write(fout);

    fout.close();

    sprintf(output_theta, "%s_minus.theta", output_name);
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

  if (state == 1 || state == 3) {
    sprintf(output_param, "%s.beta", output_name);
    fout.open(output_param);
    assert(fout.is_open());

    fout.precision(10);
    fout.unsetf(std::ios::floatfield);

    fout<< M<< std::endl;
    for (int i = 1; i <= M; ++i) transcripts[i]->write(fout);

    fout.close();

    sprintf(output_theta, "%s_plus.theta", output_name);
    fout.open(output_theta);
    assert(fout.is_open());

    fout.precision(10);
    fout.unsetf(std::ios::floatfield);

    fout<< M<< std::endl;
    fout<< prob_noise[1][0];
    for (int i = 1; i <= M; ++i) fout<< '\t'<< prob_noise[1][1] * theta[i];
    fout<< std::endl;
    
    fout.close();

    sprintf(output_rate, "%s.rate", output_name);
    sprintf(output_param, "%s.freq", output_name);
    std::ofstream fc(output_rate);
    fout.open(output_param);
    assert(fc.is_open() && fout.is_open());

    fc.precision(10);
    fc.unsetf(std::ios::floatfield);
    fout.precision(10);
    fout.unsetf(std::ios::floatfield);
    
    fc<< M<< std::endl;
    fout<< M<< std::endl;

    for (int i = 1; i <= M; ++i) {
      transcripts[i]->writeFreq(fc, fout);
      if (i < M) fc<< '\t';
      else fc<< std::endl;
    }
    
    fc.close();
    fout.close();
  }

  // write out expression results
  writeExprRes(state, output_name);

  if (verbose) printf("DMSWholeModel::write is finished!\n");
}

void DMSWholeModel::startSimulation() {
  cdf = new double[M + 1];
  cdf[0] = theta[0];
  for (int i = 1; i <= M; ++i) {
    cdf[i] = 0.0;
    if (theta[i] > 0.0) {
      assert(transcripts[i]->getEffLen() > 0);
      transcripts[i]->startSimulation();
      cdf[i] = theta[i] * transcripts[i]->getProbPass(); 
    }
    cdf[i] += cdf[i - 1];
  }
}

void DMSWholeModel::finishSimulation() {
  for (int i = 1; i <= M; ++i) 
    if (theta[i] > 0.0) transcripts[i]->finishSimulation();
  delete[] cdf;
}

void DMSWholeModel::allocateTranscriptsToThreads(int channel) {
  int id;
  MyHeap my_heap;
  std::vector<int> max_lens; // record maximum len in each thread

  // Allocate transcripts for updating
  my_heap.init(num_threads);
  std::vector<Params*> &paramsVecU = paramsVecUp[channel];
  paramsVecU.assign(num_threads, NULL);
  for (int i = 0; i < num_threads; ++i) paramsVecU[i] = new Params(i, this);

  for (int i = 1; i <= M; ++i) {
    HIT_INT_TYPE numAlign = transcripts[i]->getNumAlignments();
    if (numAlign > 0) {
      id = my_heap.getTop();
      paramsVecU[id]->trans.push_back(transcripts[i]);
      my_heap.updateTop(numAlign);
    }
  }

  // set number of transcripts per threads and delete extra paramsVec
  for (id = 0; id < num_threads; ++id) {
    paramsVecU[id]->num_trans = my_heap.getNum(id);
    if (paramsVecU[id]->num_trans == 0) break;
  }
  assert(id > 0);
  if (id < num_threads) {
    for (int i = id; i < num_threads; ++i) { 
      delete paramsVecU[i];
      paramsVecU[i] = NULL; 
    }
    paramsVecU.resize(id, NULL);
  }

  // Allocate transcripts for EM  
  my_heap.init(num_threads);
  std::vector<Params*> &paramsVecE = paramsVecEM[channel];
  paramsVecE.assign(num_threads, NULL);
  for (int i = 0; i < num_threads; ++i) paramsVecE[i] = new Params(i, this);
  max_lens.assign(num_threads, 0);

  // allocate transcripts
  for (int i = 1; i <= M; ++i) 
    if (transcripts[i]->getNumAlignments() > 0) {
      id = my_heap.getTop();
      paramsVecE[id]->trans.push_back(transcripts[i]);
      if (max_lens[id] < transcripts[i]->getLen()) max_lens[id] = transcripts[i]->getLen();
      my_heap.updateTop(transcripts[i]->getLen());
    }

  // set number of transcripts per threads and delete extra paramsVec
  for (id = 0; id < num_threads; ++id) {
    paramsVecE[id]->num_trans = my_heap.getNum(id);
    if (paramsVecE[id]->num_trans == 0) break;
  }
  assert(id > 0);
  if (id < num_threads) {
    for (int i = id; i < num_threads; ++i) { 
      delete paramsVecE[i];
      paramsVecE[i] = NULL; 
    }
    paramsVecE.resize(id, NULL);
  }

  // allocate start2 and end2 for each thread
  for (int i = 0; i < num_threads; ++i) {
    paramsVecE[i]->start2 = new double[max_lens[i] + 1];
    paramsVecE[i]->end2 = new double[max_lens[i] + 1];
    for (int j = 0; j < paramsVecE[i]->num_trans; ++j)
      paramsVecE[i]->trans[j]->setStart2andEnd2(paramsVecE[i]->start2, paramsVecE[i]->end2);
  }
}
