#include<cassert>
#include<string>
#include<fstream>

#include "Markov.hpp"
#include "Profile.hpp"
#include "QProfile.hpp"
#include "SequencingModel.hpp"

SequencingModel::SequencingModel(bool hasQual, int maxL) {
  this->hasQual = hasQual;
  markov = NULL;
  profile = NULL;
  qprofile = NULL;
  
  markov = new Markov();
  if (hasQual) qprofile = new QProfile();
  else profile = new Profile(maxL);
}

SequencingModel::~SequencingModel() {
  delete markov;
  if (hasQual) delete qprofile;
  else delete profile;
}

void SequencingModel::init() {
  markov->init();
  if (hasQual) qprofile->init();
  else profile->init();
}

void SequencingModel::collect(const SequencingModel* o) {
  markov->collect(o->markov);
  if (hasQual) qprofile->collect(o->qprofile);
  else profile->collect(o->profile);
}

void SequencingModel::finish() {
  markov->finish();
  if (hasQual) qprofile->finish();
  else profile->finish();
}

void SequencingModel::read(std::ifstream& fin) {
  std::string line;
  while (getline(fin, line)) {
    if (line.substr(0, 16) == "#SequencingModel") break;
  }
  assert(fin.good());

  markov->read(fin);
  if (hasQual) qprofile->read(fin);
  else profile->read(fin);
}

void SequencingModel::write(std::ofstream& fout) {
  fout<< "#SequencingModel"<< std::endl<< std::endl;
  markov->write(fout);
  if (hasQual) qprofile->write(fout);
  else profile->write(fout);
}

void SequencingModel::startSimulation() {
  markov->startSimulation();
  if (hasQual) qprofile->startSimulation();
  else profile->startSimulation();
}

void SequencingModel::finishSimulation() {
  markov->finishSimulation();
  if (hasQual) qprofile->finishSimulation();
  else profile->finishSimulation();
}
