#include<cmath>
#include<ctime>
#include<cstdio>
#include<cstring>
#include<cstdlib>
#include<cassert>
#include<vector>
#include<fstream>

#include "sampling.hpp"
#include "DMSTransModel.hpp"

using namespace std;

DMSTransModel *model;
Sampler *sampler;

double tot_c;
char inF[1005], outF[1005];

FILE *fi, *fo;
ifstream fin;
ofstream fout;

int primer_length, min_frag_len, max_frag_len;
double gamma_init, beta_init;
int round_minus, round_plus;

double calcL1(int len, double *a, double *b) {
  double result = 0.0;
  for (int i = 1; i <= len; ++i) result += fabs(a[i] - b[i]);
  return result;
}

double calcMean(int len, double *a) {
  double result = 0.0;
  for (int i = 1; i <= len; ++i) result += a[i];
  result /= len;
  return result;
}

int main(int argc, char* argv[]) {
  if (argc < 6 || argc > 9) {
    printf("Usage: dms_single_transcript config_file learning_from_real_data minus_channel_count.txt plus_channel_count.txt output_name\n");
    printf("       dms_single_transcript config_file simulation ground_truth.gamma ground_truth.beta output_name number_of_fragments_simulated [seed]\n");
    printf("       dms_single_transcript config_file learning_from_simulated_data minus_channel_simulated_data.dat plus_channel_simulated_data.dat output_name ground_truth_gamma.txt ground_truth_beta.txt reads_are_single_end(0 or 1)\n\n");
    printf("   config_file provides the following parameters: primer length, the minimum fragment length, maximum fragment length, mean gamma used for initializing EM, mean beta used for initialize EM, rounds of EM ran for minus channel and rounds of EM ran for plus channel. Each parameter in a separate line.\n");
    printf("   After the config_file is an action key word, it can be one of \"learing_from_real_data\", \"simulation\" or \"learning from_simulated_data\"\n\n");
    printf("   If the key word is \"learning_from_real_data\", this program try to learn parameters from a real data set, which is assumed to be single-end reads only.\n");
    printf("   The next two arguments are count vectors extracted from minus channel and plus channel. Each file contains only one file. The first value in the line is the transcript length, trans_len. Then trans_len double-precision values follows. Each value represents the number of expected read count at that transcript coordinate (from 5', 0-based). The count vectors can be extracted by script \"extract_count_vector\".\n");
    printf("   The last argument, output_name, gives the prefix for all output files. output_name.gamma, output_name.beta and output_name.theta will be produced, which represent the estimated gamma values, beta values, and theta values (c rate included). Each file contains only one line. For the first two files, the first value is the total number of parameters estimated, len (= trans_len - primer_length). Then len values follows, which are the estimated parameters (ordered from the 5' end). For the last file, the first value is the estimated rate, c.\n\n");
    printf("   If the key word is \"simulation\", this program will simulate fragment counts based on learned gamma and beta values.\n");
    printf("   ground_truth.gamma and ground_truth.beta are the parameters learned from a real data set.\n");
    printf("   output_name is the output prefix. output_name_minus.dat and output_name_plus.data will be produced, which contain the simulated fragments. The format of these files is as follows: Each line contains a simulated fragment. In each line, the first value is this fragment's start position (0-based, from 5' end) and the second value is the fragment length.\n");
    printf("   number_of_fragments_simulated sets the number of fragments we want to simulate.\n");
    printf("   [seed] is an optional argument. It sets the seed used for simulation. A same seed will always result in a same simulation.\n\n");
    printf("   If the key word is \"learning_from_simulated_data\", this program will learn parameters from simulated data sets.\n");
    printf("   minus_channel_simulated_data.dat and plus_channel_simulated_data.data give the simulated data themselves.\n");
    printf("   output_name is defined above. output_name.gamma, output_name.beta and output_name.theta will be generated.\n");
    printf("   ground_truth_gamma.txt and ground_truth_beta.txt provides the ground truth parameters, which should be learned from real data sets.\n");
    printf("   reads_are_single_end(0 or 1) tells if we should view the simulated fragments as single-end reads or paired-end reads. 0 means reads are paired-end and 1 means reads are single-end.\n\n");
    exit(-1);
  }

  fi = fopen(argv[1], "r");
  assert(fscanf(fi, "%d %d %d %lf %lf %d %d", &primer_length, &min_frag_len, &max_frag_len, &gamma_init, &beta_init, &round_minus, &round_plus) == 7);
  DMSTransModel::setGlobalParams(primer_length, min_frag_len, max_frag_len, gamma_init, beta_init);

  sampler = NULL;
  model = NULL;

  if (!strcmp(argv[2], "learning_from_real_data")) {
    int trans_len;
    vector<double> counts;

    fi = fopen(argv[3], "r");
    fscanf(fi, "%d", &trans_len);
    counts.assign(trans_len, 0);
    tot_c = 0.0;
    for (int i = 0; i < trans_len; ++i) {
      fscanf(fi, "%lf", &counts[i]);
      tot_c += counts[i];
    }

    printf("MINUS_TOTAL_COUNTS = %.2f\n", tot_c);

    model = new DMSTransModel(true, trans_len);

    model->init();
    for (int i = 0; i < trans_len; ++i) 
      if (counts[i] > 0.0) model->update(i, counts[i]);
    fclose(fi);

    model->calcAuxiliaryArrays();
    model->EM(tot_c, round_minus);

    printf("MINUS_PROB_PASS = %.10g, Loglik = %.10g\n", model->getProbPass(), model->calcLogLik());

    sprintf(outF, "%s.gamma", argv[5]);
    fout.open(outF);
    model->write(fout);
    fout.close();

    sprintf(inF, "%s.gamma", argv[5]);
    fin.open(inF);
    model->read(fin);
    fin.close();

    fi = fopen(argv[4], "r");
    fscanf(fi, "%d", &trans_len);
    tot_c = 0.0;
    for (int i = 0; i < trans_len; ++i) {
      fscanf(fi, "%lf", &counts[i]);
      tot_c += counts[i];
    }

    printf("PLUS_TOT_C = %.2f\n", tot_c);

    model->init();
    for (int i = 0; i < trans_len; ++i) 
      if (counts[i] > 0.0) model->update(i, counts[i]);
    fclose(fi);

    model->calcAuxiliaryArrays();
    model->EM(tot_c, round_plus);

    printf("PLUS_PROB_PASS = %.10g, Loglik = %.10g\n", model->getProbPass(), model->calcLogLik());

    sprintf(outF, "%s.beta", argv[5]);
    fout.open(outF);
    model->write(fout);
    fout.close();

    sprintf(outF, "%s.theta", argv[5]);
    fout.open(outF);
    model->writeTheta(fout);
    fout.close();

    delete model;
  }
  else if (!strcmp(argv[2], "learning_from_simulated_data")) {
    int len;
    double *gt_minus = NULL, *gt_plus = NULL;
    //    double *old_minus = NULL, *old_plus = NULL;
    bool isSE = (atoi(argv[8]) == 1);

    //    int best_round = 0;
    //    double best_l1 = 1e100, l1_value = 1e100;

    // Load ground truth
    fi = fopen(argv[6], "r");
    fscanf(fi, "%d", &len);
    gt_minus = new double[len + 1];
    gt_minus[0] = 0.0;
    for (int i = 1; i <= len; ++i) fscanf(fi, "%lf", &gt_minus[i]);
    fclose(fi);

    fi = fopen(argv[7], "r");
    fscanf(fi, "%d", &len);
    gt_plus = new double[len + 1];
    gt_plus[0] = 0.0;
    for (int i = 1; i <= len; ++i) fscanf(fi, "%lf", &gt_plus[i]);
    fclose(fi);

    int pos, fragment_length;
    model = new DMSTransModel(true, len + DMSTransModel::get_primer_length());
    model->init();

    fi = fopen(argv[3], "r");
    tot_c = 0.0;
    while (fscanf(fi, "%d %d", &pos, &fragment_length) == 2) {
      if (isSE) model->update(pos, 1.0);
      else model->update(pos, fragment_length, 1.0);      
      ++tot_c;
    }
    fclose(fi);

    model->calcAuxiliaryArrays();
    model->EM(tot_c, round_minus);
    printf("L1_diff = %.10g, Loglik = %.10g, Probability_of_a_fragment_passing_the_size_selection_step =  %.10g\n", calcL1(len, gt_minus, model->getGamma()), model->calcLogLik(), model->getProbPass());
    /*
    best_round = 0;
    best_l1 = 1e100;
    
    old_minus = new double[len + 1];
    memcpy(old_minus, model->getGamma(), sizeof(double) * (len + 1));
    for (int round = 1; round <= 2000; ++round) {
      model->EM2(tot_c, 1);
      l1_value = calcL1(len, gt_minus, model->getGamma());
      if (best_l1 > l1_value) { best_l1 = l1_value; best_round = round; }
      printf("ROUND = %d, L1_diff = %.10g, Diff = %.10g\n", round, l1_value, calcL1(len, old_minus, model->getGamma()));
      memcpy(old_minus, model->getGamma(), sizeof(double) * (len + 1));
     }
  
    printf("Best round = %d, best_l1 = %.10g\n", best_round, best_l1);
    delete[] old_minus;
    */

    sprintf(outF, "%s.gamma", argv[5]);
    fout.open(outF);
    model->write(fout);
    fout.close();

    sprintf(inF, "%s.gamma", argv[5]);
    fin.open(inF);
    model->read(fin);
    fin.close();

    fi = fopen(argv[4], "r");
    tot_c = 0.0;
    model->init();
    while (fscanf(fi, "%d %d", &pos, &fragment_length) == 2) {
      if (isSE) model->update(pos, 1.0);
      else model->update(pos, fragment_length, 1.0);
      ++tot_c;
    }
    fclose(fi);

    model->calcAuxiliaryArrays();
    model->EM(tot_c, round_plus);
    printf("L1_diff = %.10g, Loglik = %.10g, Prob_pass = %.10g\n", calcL1(len, gt_plus, model->getBeta()), model->calcLogLik(), model->getProbPass());
    /*
    best_round = 0;
    best_l1 = 1e100;

    old_plus = new double[len + 1];
    memcpy(old_plus, model->getBeta(), sizeof(double) * (len + 1));
    for (int round = 1; round <= 2000; ++round) {
      model->EM(tot_c, 1);
      l1_value = calcL1(len, gt_plus, model->getBeta());
      if (best_l1 > l1_value) { best_l1 = l1_value; best_round = round; }      
      printf("ROUND = %d, L1_diff = %.10g, DIFF = %.10g\n", round, l1_value, calcL1(len, old_plus, model->getBeta()));
      memcpy(old_plus, model->getBeta(), sizeof(double) * (len + 1));
    }
    delete[] old_plus;
    printf("Best_round = %d, Best_l1 = %.10g\n", best_round, best_l1);
    */

    sprintf(outF, "%s.beta", argv[5]);
    fout.open(outF);
    model->write(fout);
    fout.close();

    sprintf(outF, "%s.theta", argv[5]);
    fout.open(outF);
    model->writeTheta(fout);
    fout.close();

    delete model;
    delete[] gt_minus;
    delete[] gt_plus;
  }
  else {
    assert(!strcmp(argv[2], "simulation"));
    int num_reads = atoi(argv[6]);

    seedType seed = 0; 
    if (argc == 8) {
      int str_len = strlen(argv[argc - 1]);
      for (int i = 0; i < str_len; ++i) seed = seed * 10 + (argv[argc - 1][i] - '0');
    }
    else seed = time(NULL);
    sampler = new Sampler(seed);
    
    model = new DMSTransModel(false);

    // Simulate minus channel
    fin.open(argv[3]);
    model->read(fin);
    fin.close();

    sprintf(outF, "%s_minus.dat", argv[5]);
    fo = fopen(outF, "w");
    int pos, fragment_length;

    for (int i = 0; i < num_reads; ++i) {
      model->simulate(sampler, pos, fragment_length);
      fprintf(fo, "%d\t%d\n", pos, fragment_length);
    }
    fclose(fo);

    printf("MINUS SIMULATION FIN!\n");

    // Simulate plus channel
    fin.open(argv[4]);
    model->read(fin);
    fin.close();
    
    sprintf(outF, "%s_plus.dat", argv[5]);
    fo = fopen(outF, "w");
    for (int i = 0; i < num_reads; ++i) {
      model->simulate(sampler, pos, fragment_length);
      fprintf(fo, "%d\t%d\n", pos, fragment_length);
    }
    fclose(fo);

    printf("PLUS SIMULATION FIN!\n");

    delete sampler;
    delete model;
  }

  return 0;
}
