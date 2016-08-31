/* Copyright (c) 2016
   Bo Li (University of California, Berkeley)
   bli25@berkeley.edu

   This program is free software; you can redistribute it and/or
   modify it under the terms of the GNU General Public License as
   published by the Free Software Foundation; either version 3 of the
   License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
   General Public License for more details.   

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
   USA
*/

#include <cassert>
#include <string>
#include <fstream>

#include "FragLenDist.hpp"
#include "SequencingModel.hpp"
#include "PROBerReadModel_iCLIP.hpp"

PROBerReadModel_iCLIP::PROBerReadModel_iCLIP() : model_type(-1), fld(NULL), seqmodel(NULL) {}

PROBerReadModel_iCLIP::PROBerReadModel_iCLIP(int model_type, int max_len) : model_type(model_type) {
  fld = NULL; seqmodel = NULL;
  if (model_type >= 2) fld = new FragLenDist(max_len);
  seqmodel = new SequencingModel((model_type & 1), max_len);
  seqmodel->init(); // initialize seqmodel for updates
}

PROBerReadModel_iCLIP::~PROBerReadModel_iCLIP() {
  if (fld != NULL) delete fld;
  if (seqmodel != NULL) delete seqmodel;
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

  // Read FragLenDist if model_type >= 2
  if (model_type >= 2) {
    fld = new FragLenDist();
    fld->read(fin);
  }

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

  // Write FragLenDist if model_type >= 2
  if (model_type >= 2) fld->write(fout);
  
  // Write sequencing model
  seqmodel->write(fout);

  fout.close();

  if (verbose) printf("PROBerReadModel::write finished!\n");
}


