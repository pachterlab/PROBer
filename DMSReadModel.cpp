#include<string>
#include<fstream>
#include<algorithm>

#include "utils.h"
#include "sampling.hpp"

#include "Refs.hpp"

#include "MateLenDist.hpp"
#include "SequencingModel.hpp"
#include "NoiseProfile.hpp"
#include "QualDist.hpp"

#include "DMSReadModel.hpp"

DMSReadModel::DMSReadModel(int model_type, Refs* refs) : model_type(model_type), refs(refs) {
  mld1 = new MateLenDist();
  if (model_type >= 2) mld2 = new MateLenDist();
  if (model_type == 1 || model_type == 3) qd = new QualDist();
  npro = new NoiseProfile(true);

  seqmodel = NULL;

  max_len = 0;
  loglik = 0.0;
  sampler = NULL;
}

DMSReadModel::DMSReadModel(DMSReadModel* master_model) {
  model_type = master_model->model_type;
  mld1 = mld2 = NULL;
  qd = NULL;
  max_len = master_model->max_len;
  loglik = 0.0;

  npro = new NoiseProfile();
  seqmodel = new SequencingModel((model_type == 1 || model_type == 3), max_len);

  refs = master_model->refs;
  sampler = NULL;
}

DMSReadModel::DMSReadModel(Refs* refs, Sampler* sampler) : refs(refs), sampler(sampler) {
  model_type = -1; // not initilaized yet

  mld1 = mld2 = NULL;
  qd = NULL;
  seqmodel = NULL;
  npro = NULL;

  max_len = 0;
  loglik = 0.0;
}

DMSReadModel::~DMSReadModel() {
  if (mld1 != NULL) delete mld1;
  if (mld2 != NULL) delete mld2;
  if (qd != NULL) delete qd;
  if (seqmodel != NULL) delete seqmodel;
  if (npro != NULL) delete npro;

  refs = NULL;
  sampler = NULL;
}

void DMSReadModel::finish_preprocess() {
  loglik = 0.0;
  mld1->finish();
  loglik += mld1->getLogP();
  if (model_type >= 2) {
    mld2->finish();
    loglik += mld2->getLogP();
  }
  if (model_type & 1) loglik += qd->finish();
  npro->calcInitParams();

  max_len = mld1->getMaxL();
  if (model_type >= 2 && max_len < mld2->getMaxL()) max_len = mld2->getMaxL();
  seqmodel = new SequencingModel((model_type == 1 || model_type == 3), max_len);
}

void DMSReadModel::init() {
  seqmodel->init();
  npro->init();
}

void DMSReadModel::collect(DMSReadModel* o) {
  seqmodel->collect(o->seqmodel);
  npro->collect(o->npro);
}

void DMSReadModel::finish() {
  seqmodel->finish();
  npro->finish();
}

void DMSReadModel::read(const char* modelF) {
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

  // Read mate length distribution 1
  mld1 = new MateLenDist();
  mld1->read(fin);

  // Read mate length distribution 2
  if (model_type >= 2) {
    mld2 = new MateLenDist();
    mld2->read(fin);
  }

  // Read QualDist
  if (model_type & 1) {
    qd = new QualDist();
    qd->read(fin);
  }

  // Read sequencing model
  seqmodel = new SequencingModel((model_type & 1));
  seqmodel->read(fin);

  // Read noise profile model
  npro = new NoiseProfile();
  npro->read(fin);

  fin.close();

  if (verbose) printf("DMSReadModel::read finished!\n");
}

void DMSReadModel::write(const char* modelF) {
  std::ofstream fout(modelF);
  assert(fout.is_open());

  fout.precision(10);
  fout.unsetf(std::ios::floatfield);

  // Write model type
  fout<< "#Model Type: 0, SE, no qual; 1, SE, qual; 2, PE, no qual; 3 PE, qual"<< std::endl;
  fout<< model_type<< std::endl<< std::endl;

  // Write mate length distribution 1
  fout<< "#Mate length distribution 1"<< std::endl;
  mld1->write(fout);

  // Write mate length distribution 2
  if (model_type >= 2) {
    fout<< "#Mate length distribution 2"<< std::endl;
    mld2->write(fout);
  }

  // Write QualDist
  if (model_type & 1) qd->write(fout);

  // Write sequencing model
  seqmodel->write(fout);

  // Write noise profile
  npro->write(fout);

  fout.close();

  if (verbose) printf("DMSReadModel::write finished!\n");
}

void DMSReadModel::startSimulation() {
  if (model_type & 1) qd->startSimulation();
  seqmodel->startSimulation();
  npro->startSimulation();
}

void DMSReadModel::finishSimulation() {
  if (model_type & 1) qd->finishSimulation();
  seqmodel->finishSimulation();
  npro->finishSimulation();
}
