#include<cassert>
#include<string>
#include<vector>
#include<fstream>
#include<stdint.h>

#include "RefSeq.hpp"

RefSeq::RefSeq() {
  fullLen = totLen = 0;
  name = ""; seq = "";
  fmasks.clear();
  dir = 0;
}

//Constructor , seq : the forward strand of the reference
//tag does not contain ">"
//polyALen : length of polyA tail we add
RefSeq::RefSeq(const std::string& name, const std::string& seq, int polyALen) {
  fullLen = seq.length();
  totLen = fullLen + polyALen;
  this->name = name;
  this->seq = seq;
  this->seq.append(polyALen, 'A');
  
  assert(fullLen > 0 && totLen >= fullLen);
  
  int len = (fullLen - 1) / NBITS + 1;
  fmasks.assign(len, 0);

  // Set mask if poly(A) tail is added
  // This part will be removed once we have a good way to add mask sequences
  if (polyALen > 0) {
    int OLEN = 25; // last 24 bases are masked as no read should align to
    for (int i = std::max(fullLen - OLEN + 1, 0); i < fullLen; i++) setMask(i);
  }

  // Default direction is '+'
  dir = '+';
}

RefSeq::RefSeq(const RefSeq& o) {
  fullLen = o.fullLen;
  totLen = o.totLen;
  name = o.name;
  seq = o.seq;
  fmasks = o.fmasks;
  dir = o.dir;
}

RefSeq& RefSeq::operator= (const RefSeq &rhs) {
  if (this != &rhs) {
    fullLen = rhs.fullLen;
    totLen = rhs.totLen;
    name = rhs.name;
    seq = rhs.seq;
    fmasks = rhs.fmasks;
    dir = rhs.dir;
  }
  
  return *this;
}

//internal read; option 0 : read all 1 : do not read seqences
bool RefSeq::read(std::ifstream& fin, int option) {
  std::string line;
  
  if (!(fin>>fullLen>>totLen)) return false;
  assert(fullLen > 0 && totLen >= fullLen);
  getline(fin, line);
  if (!getline(fin, name)) return false;
  if (!getline(fin, seq)) return false;
  
  int len = (fullLen - 1) / NBITS + 1; // assume each cell contains NBITS bits
  fmasks.assign(len, 0);
  for (int i = 0; i < len; i++)
    if (!(fin>>fmasks[i])) return false;
  getline(fin, line);
  
  assert(option == 0 || option == 1);
  if (option == 1) { seq = ""; }

  dir = '+'; // Default is '+'

  return true;
}

//write to file in "internal" format
void RefSeq::write(std::ofstream& fout) {
  fout<<fullLen<<" "<<totLen<<std::endl;
  fout<<name<<std::endl;
  fout<<seq<<std::endl;
  
  int len = fmasks.size();
  for (int i = 0; i < len - 1; i++) fout<<fmasks[i]<<" ";
  fout<<fmasks[len - 1]<<std::endl;
}

const std::vector<uint32_t> RefSeq::mask_codes = RefSeq::init_mask_code();

std::vector<uint32_t> RefSeq::init_mask_code() {
  std::vector<uint32_t> vec(NBITS);
  for (int i = 0; i < NBITS; i++) vec[i] = 1 << i;
  return vec;
}

const std::vector<int> RefSeq::base2code = init_base2code();

std::vector<int> RefSeq::init_base2code() {
  std::vector<int> vec(128, -1);
  vec['a'] = vec['A'] = 0;
  vec['c'] = vec['C'] = 1;
  vec['g'] = vec['G'] = 2;
  vec['t'] = vec['T'] = 3;
  vec['n'] = vec['N'] = 4;
  
  return vec;
}

const std::vector<int> RefSeq::rbase2code = RefSeq::init_rbase2code();

std::vector<int> RefSeq::init_rbase2code() {
  std::vector<int> vec(128, -1);
  vec['a'] = vec['A'] = 3;
  vec['c'] = vec['C'] = 2;
  vec['g'] = vec['G'] = 1;
  vec['t'] = vec['T'] = 0;
  vec['n'] = vec['N'] = 4;
  
  return vec;
}
