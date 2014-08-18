#include<cassert>
#include<string>
#include<fstream>

#include "RSPD.hpp"

void RSPD::read(std::ifstream& fin) {
  std::string line;
  while (getline(fin, line)) {
    if (line.substr(0, 5) == "#RSPD") break;
  }
  assert((fin>> estRSPD) && !estRSPD);
  getline(fin, line);
}

void RSPD::write(std::ofstream& fout) {
  fout<< "#RSPD"<< std::endl;
  fout<< estRSPD<< std::endl<< std::endl;
}
