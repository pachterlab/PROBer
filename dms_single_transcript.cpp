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

double calcL1(int len, double *a, double *b) {
  double result = 0.0;
  for (int i = 1; i <= len; ++i) result += fabs(a[i] - b[i]);
  return result;
}

int main(int argc, char* argv[]) {
  if (argc < 5 || argc > 8) {
    printf("Usage: dms_single_transcript learning_from_real_data minus_channel_count.txt plus_channel_count.txt output_name\n");
    printf("       dms_single_transcript simulation gamma.txt beta.txt output_name number_of_reads_simulated [seed]\n");
    printf("       dms_single_transcript learning_from_simulated_data minus_channel_simulated_data.dat plus_channel_simulated_data.dat output_name ground_truth_gamma.txt ground_truth_beta.txt reads_are_single_end(0 or 1)\n");
    exit(-1);
  }

  //DMSTransModel::setGlobalParams(6, 25, 45);
  DMSTransModel::setGlobalParams(6, 150, 650);

  model = NULL;
  sampler = NULL;

  if (!strcmp(argv[1], "learning_from_real_data")) {
    int trans_len;
    vector<double> counts;

    fi = fopen(argv[2], "r");
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
    model->EM(tot_c, 200);

    printf("MINUS_PROB_PASS = %.10g\n", model->getProbPass());

    sprintf(outF, "%s.gamma", argv[4]);
    fout.open(outF);
    model->write(fout);
    fout.close();

    sprintf(inF, "%s.gamma", argv[4]);
    fin.open(inF);
    model->read(fin);
    fin.close();

    fi = fopen(argv[3], "r");
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
    model->EM(tot_c, 200);

    printf("PLUS_PROB_PASS = %.10g\n", model->getProbPass());

    sprintf(outF, "%s.beta", argv[4]);
    fout.open(outF);
    model->write(fout);
    fout.close();

    delete model;
  }
  else if (!strcmp(argv[1], "learning_from_simulated_data")) {
    int len;
    double *gt_minus = NULL, *gt_plus = NULL;
    double *old_minus = NULL, *old_plus = NULL;
    bool isSE = (atoi(argv[7]) == 1);

    int best_round = 0;
    double best_l1 = 1e100, l1_value = 1e100;

    // Load ground truth
    fi = fopen(argv[5], "r");
    fscanf(fi, "%d", &len);
    gt_minus = new double[len + 1];
    gt_minus[0] = 0.0;
    for (int i = 1; i <= len; ++i) fscanf(fi, "%lf", &gt_minus[i]);
    fclose(fi);

    fi = fopen(argv[6], "r");
    fscanf(fi, "%d", &len);
    gt_plus = new double[len + 1];
    gt_plus[0] = 0.0;
    for (int i = 1; i <= len; ++i) fscanf(fi, "%lf", &gt_plus[i]);
    fclose(fi);

    int pos, fragment_length;
    model = new DMSTransModel(true, len + DMSTransModel::get_primer_length());
    model->init();

    fi = fopen(argv[2], "r");
    tot_c = 0.0;
    while (fscanf(fi, "%d %d", &pos, &fragment_length) == 2) {
      if (isSE) model->update(pos, 1.0);
      else model->update(pos, fragment_length, 1.0);      
      ++tot_c;
    }
    fclose(fi);

    //    model->calcAuxiliaryArrays();
    //model->EM(tot_c, 420);
    //printf("%.10g %.10g\n", calcL1(len, gt_minus, model->getGamma()), model->getProbPass());

     
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
    exit(-1);

    sprintf(outF, "%s.gamma", argv[4]);
    fout.open(outF);
    model->write(fout);
    fout.close();

    sprintf(inF, "%s.gamma", argv[4]);
    fin.open(inF);
    model->read(fin);
    fin.close();

    fi = fopen(argv[3], "r");
    tot_c = 0.0;
    model->init();
    while (fscanf(fi, "%d %d", &pos, &fragment_length) == 2) {
      if (isSE) model->update(pos, 1.0);
      else model->update(pos, fragment_length, 1.0);
      ++tot_c;
    }
    fclose(fi);

    model->calcAuxiliaryArrays();
    //model->EM(tot_c, 1041);
    //printf("%.10g %.10g\n", calcL1(len, gt_plus, model->getBeta()), model->getProbPass());
    
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

    sprintf(outF, "%s.beta", argv[4]);
    fout.open(outF);
    model->write(fout);
    fout.close();

    delete model;
    delete[] gt_minus;
    delete[] gt_plus;
  }
  else {
    assert(!strcmp(argv[1], "simulation"));
    int num_reads = atoi(argv[5]);

    seedType seed = 0;
    if (argc == 7) {
      int str_len = strlen(argv[6]);
      for (int i = 0; i < str_len; ++i) seed = seed * 10 + (argv[6][i] - '0');
    }
    else seed = time(NULL);
    
    sampler = new Sampler(seed);

    model = new DMSTransModel(false);

    // Simulate minus channel
    fin.open(argv[2]);
    model->read(fin);
    fin.close();

    sprintf(outF, "%s_minus.dat", argv[4]);
    fo = fopen(outF, "w");
    int pos, fragment_length;

    for (int i = 0; i < num_reads; ++i) {
      model->simulate(sampler, pos, fragment_length);
      fprintf(fo, "%d\t%d\n", pos, fragment_length);
    }
    fclose(fo);

    printf("MINUS SIMULATION FIN!\n");

    // Simulate plus channel
    fin.open(argv[3]);
    model->read(fin);
    fin.close();
    
    sprintf(outF, "%s_plus.dat", argv[4]);
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
