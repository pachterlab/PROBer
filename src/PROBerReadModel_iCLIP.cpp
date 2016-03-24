#include<cassert>
#include<string>
#include<fstream>

#include "utils.h"
#include "SequencingModel.hpp"
#include "PROBerReadModel_iCLIP.hpp"

PROBerReadModel_iCLIP::PROBerReadModel_iCLIP(int model_type, int max_len) : model_type(model_type) {
  seqmodel = new SequencingModel((model_type & 1), max_len);
}

PROBerReadModel_iCLIP::~PROBerReadModel_iCLIP() {
  delete seqmodel;
}

void PROBerReadModel_iCLIP::read(const char* modelF) {
  std::string line;
  std::ifstream fin(modelF);
  assert(fin.is_open());

  // Read model type
  while (getline(fin, line)) {
    if (line.substr(0, 11) == "#Model Type") break;
  }
  assert(fin.good());
  assert(fin>> model_type);
  getline(fin, line);

  // Read sequencing model
  seqmodel = new SequencingModel((model_type & 1));
  seqmodel->read(fin);

  fin.close();

  if (verbose) printf("PROBerReadModel_iCLIP::read finished!\n");
}

void PROBerReadModel_iCLIP::write(const char* modelF) {
  std::ofstream fout(modelF);
  assert(fout.is_open());

  fout.precision(10);
  fout.unsetf(std::ios::floatfield);

  // Write model type
  fout<< "#Model Type: 0, SE, no qual; 1, SE, qual; 2, PE, no qual; 3 PE, qual"<< std::endl;
  fout<< model_type<< std::endl<< std::endl;

  // Write sequencing model
  seqmodel->write(fout);

  fout.close();

  if (verbose) printf("PROBerReadModel::write finished!\n");
}


