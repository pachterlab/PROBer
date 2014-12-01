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

DMSWholeModel::DMSWholeModel(const char* config_file, const Transcripts* trans, int num_threads, int read_length) {
  // set DMSTransModel static member values
  int primer_length, min_frag_len, max_frag_len;
  double gamma_init, beta_init;

  FILE *fi = fopen(config_file, "r");
  assert(fi != NULL);
  assert(fscanf(fi, "%d %d %d %lf %lf", &primer_length, &min_frag_len, &max_frag_len, &gamma_init, &beta_init) == 5);
  fclose(fi);

  DMSTransModel::setGlobalParams(primer_length, min_frag_len, max_frag_len, gamma_init, beta_init, read_length);

  readGamma = true;
  this->num_threads = 0;

  M = 0; 
  theta.clear();
  transcripts.clear();

  N_tot = 0.0;
  prob_pass = 0.0;
  counts.clear();
  unobserved.clear();

  threads.clear();
  paramsVec.clear();
  
  if (trans != NULL) {
    assert(num_threads >= 1);
    this->num_threads = num_threads;

    M = trans->getM();
    theta.assign(M + 1, 0.0);
    transcripts.assign(M + 1, NULL);
    for (int i = 1; i <= M; ++i) {
      const Transcript& tran = trans->getTranscriptAt(i);
      transcripts[i] = new DMSTransModel(true, tran.getTranscriptID(), tran.getLength());
    }
    
    counts.assign(M + 1, 0.0);
    unobserved.assign(M + 1, 0.0);
  }

  cdf = NULL;
}

DMSWholeModel::~DMSWholeModel() {
  assert(transcripts[0] == NULL);
  for (int i = 1; i <= M; ++i) delete transcripts[i];

  if (paramsVec.size() > 0) {
    pthread_attr_destroy(&attr);
    for (int i = 0; i < (int)paramsVec.size(); ++i) delete paramsVec[i];
  }
}

void DMSWholeModel::init_for_EM() {
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
  threads.assign(num_threads, pthread_t());

  allocateTranscriptsToThreads();

  // calcAuxiliaryArrays
  // create threads
  for (int i = 0; i < num_threads; ++i) {
    rc = pthread_create(&threads[i], &attr, run_calcAuxiliaryArrays_per_thread, (void*)paramsVec[i]);
    pthread_assert(rc, "pthread_create", "Cannot create thread " + itos(i) + " (numbered from 0) for run_calcAuxiliaryArrays_per_thread!");
  }
  // join threads
  for (int i = 0; i < num_threads; ++i) {
    rc = pthread_join(threads[i], NULL);
    pthread_assert(rc, "pthread_join", "Cannot join thread " + itos(i) + " (numbered from 0) for run_calcAuxiliaryArrays_per_thread!");
  }  

  // Initialize theta
  double denom = 1.0; // transcript "0" will always count
  for (int i = 1; i <= M; ++i) 
    if (!isZero(counts[i])) ++denom;
  for (int i = 0; i <= M; ++i)
    theta[i] = ((i == 0 || !isZero(counts[i])) ? 1.0 / denom : 0.0);

  calcProbPass();
}

void DMSWholeModel::update(double count0) {
  // create threads
  for (int i = 0; i < num_threads; ++i) {
    rc = pthread_create(&threads[i], &attr, run_makeUpdates_per_thread, (void*)(paramsVec[i]));
    pthread_assert(rc, "pthread_create", "Cannot create thread " + itos(i) + " (numbered from 0) for run_makeUpdates_per_thread!");
  }
  // join threads
  for (int i = 0; i < num_threads; ++i) {
    rc = pthread_join(threads[i], NULL);
    pthread_assert(rc, "pthread_join", "Cannot join thread " + itos(i) + " (numbered from 0) for run_makeUpdates_per_thread!");
  }    

  // Set counts
  counts[0] = count0;
  for (int i = 1; i <= M; ++i)
    counts[i] = transcripts[i]->getNobs();
}

void DMSWholeModel::runEM(double count0) {
  // Run EM
  double N_obs, sum;

  // Update counts
  update(count0);

  // Calculate total number of observed reads
  N_obs = 0.0;
  for (int i = 0; i <= M; ++i) 
    if (isZero(counts[i])) counts[i] = 0.0;
    else {
      N_obs += counts[i];
    }
  N_tot = N_obs / prob_pass;

  // Calculate expected hidden reads to each transcript
  for (int i = 0; i <= M; ++i) {
    unobserved[i] = 0.0;
    // If no counts, force the unobserved counts to be 0!
    if (i > 0 && !isZero(counts[i])) unobserved[i] = N_tot * theta[i] * (1.0 - transcripts[i]->getProbPass()); 
  }
  
  // Estimate new gamma/beta parameters
  // create threads
  for (int i = 0; i < num_threads; ++i) {
    rc = pthread_create(&threads[i], &attr, run_EM_step_per_thread, (void*)(paramsVec[i]));
    pthread_assert(rc, "pthread_create", "Cannot create thread " + itos(i) + " (numbered from 0) for run_EM_step_per_thread!");
  }
  // join threads
  for (int i = 0; i < num_threads; ++i) {
    rc = pthread_join(threads[i], NULL);
    pthread_assert(rc, "pthread_join", "Cannot join thread " + itos(i) + "(numbered from 0) for run_EM_step_per_thread!");
  }

  // Estimate new theta
  sum = 0.0;
  for (int i = 0; i <= M; ++i) {
    theta[i] = counts[i] + unobserved[i];
    sum += theta[i];
  }
  assert(!isZero(sum));
  for (int i = 0; i <= M; ++i) theta[i] /= sum;

  calcProbPass();
  N_tot = N_obs / prob_pass;

  if (verbose) printf("DMSWholeModel EM: prob_pass = %.10g\n", prob_pass);
}

void DMSWholeModel::read(const char* input_name) {
  char input_param[STRLEN];
  char input_theta[STRLEN];
  std::ifstream fin;

  int tmp_M;

  if (readGamma) {
    sprintf(input_param, "%s.gamma", input_name);
    fin.open(input_param);
    assert(fin.is_open());

    assert(fin>> tmp_M);
    if (M == 0) {
      M = tmp_M;
      theta.assign(M + 1, 0.0);
      transcripts.assign(M + 1, NULL);
      for (int i = 1; i <= M; ++i) transcripts[i] = new DMSTransModel(false);
      tmp_M = 0; // To mark that this model is not used for learning
    }
    else assert(M == tmp_M);

    for (int i = 1; i <= M; ++i) transcripts[i]->read(fin);

    fin.close();

    if (tmp_M == 0) {
      // Loading expression of minus chanel for the sake of simulation 
      sprintf(input_theta, "%s_minus.theta", input_name);
      fin.open(input_theta);
      assert(fin.is_open());
      
      assert((fin>> tmp_M) && (tmp_M == M));
      for (int i = 0; i <= M; ++i) assert(fin>> theta[i]);
      fin.close();
    }

    readGamma = false;
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

void DMSWholeModel::setDefault() {
  for (int i = 1; i <= M; ++i) transcripts[i]->setDefault();
}

void DMSWholeModel::writeExprRes(const char* exprF) {
  std::ofstream fout(exprF);
  assert(fout.is_open());

  fout.precision(2);
  fout.setf(std::ios::fixed, std::ios::floatfield);

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

  fout<< "transcript_id\tlength\teffective_length\texpected_count\tTPM\tFPKM"<< std::endl;
  for (int i = 1; i <= M; ++i) {
    fout<< transcripts[i]->getName()<< '\t'<< transcripts[i]->getLen() + transcripts[i]->get_primer_length()<< '\t'<< transcripts[i]->getLen() + 1<< '\t'<< counts[i]<< '\t'<< tpm[i]<< '\t'<< fpkm[i]<< '\t'<< std::endl;
  }

  fout.close();
}

void DMSWholeModel::write(const char* output_name) {
  char output_param[STRLEN];
  char output_theta[STRLEN];
  char output_rate[STRLEN];
  char exprF[STRLEN];

  std::ofstream fout;

  if (readGamma) {
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
    for (int i = 0; i < M; ++i) fout<< theta[i]<< '\t';
    fout<< theta[M]<< std::endl;
    
    fout.close();

    sprintf(exprF, "%s_minus.expr", output_name);
    writeExprRes(exprF);
  }
  else {
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
    for (int i = 0; i < M; ++i) fout<< theta[i]<< '\t';
    fout<< theta[M]<< std::endl;
    
    fout.close();

    sprintf(exprF, "%s_plus.expr", output_name);
    writeExprRes(exprF);

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

void DMSWholeModel::allocateTranscriptsToThreads() {
  int id;
  MyHeap my_heap;
  std::vector<int> max_lens; // record maximum len in each thread

  my_heap.init(num_threads);
  paramsVec.assign(num_threads, NULL);
  for (int i = 0; i < num_threads; ++i) paramsVec[i] = new Params(i, this);
  max_lens.assign(num_threads, 0);

  // allocate transcripts
  for (int i = 1; i <= M; ++i) 
    if (!isZero(counts[i])) {
      id = my_heap.getTop();
      ++paramsVec[id]->num_trans;
      paramsVec[id]->trans.push_back(transcripts[i]);
      paramsVec[id]->origin_ids.push_back(i);
      if (max_lens[id] < transcripts[i]->getLen()) max_lens[id] = transcripts[i]->getLen();
      my_heap.updateTop(transcripts[i]->getLen());
    }

  // trim empty threads  
  while (num_threads > 0 && paramsVec[num_threads - 1]->num_trans == 0) {
    --num_threads;
    delete paramsVec[num_threads];
    paramsVec[num_threads] = NULL;
  }
  paramsVec.resize(num_threads, NULL);  

  assert(num_threads > 0);

  // allocate start2 and end2 for each thread
  for (int i = 0; i < num_threads; ++i) {
    paramsVec[i]->start2 = new double[max_lens[i] + 1];
    paramsVec[i]->end2 = new double[max_lens[i] + 1];
    for (int j = 0; j < paramsVec[i]->num_trans; ++j)
      paramsVec[i]->trans[j]->setStart2andEnd2(paramsVec[i]->start2, paramsVec[i]->end2);
  }
}
