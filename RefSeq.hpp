#ifndef REFSEQ_H_
#define REFSEQ_H_

#include<cassert>
#include<fstream>
#include<string>
#include<vector>
#include<stdint.h>

//Each Object can only be used once
class RefSeq {
public:
  RefSeq();
  RefSeq(const std::string&, const std::string&, int);
  RefSeq(const RefSeq& o);
  RefSeq& operator= (const RefSeq&);
  
  bool read(std::ifstream&, int  = 0);
  void write(std::ofstream&);

  int getFullLen() const { return fullLen; }  

  int getTotLen() const { return totLen; }

  const std::string& getName() const { return name; }

  std::string getSeq() const { return seq; }

  void setDir(char dir) { 
    assert(dir == '+' || dir == '-');
    this->dir = dir; 
  }

  char baseAt(int pos) const {
    assert(pos >= 0 && pos < totLen);
    return (dir == '+' ? seq[pos] : seq[totLen - pos - 1]);
  }

  int baseCodeAt(int pos) const {
    assert(pos >= 0 && pos < totLen);
    return (dir == '+' ? base2code[seq[pos]] : rbase2code[seq[totLen - pos - 1]]);
  }
  
  bool getMask(int pos) const {
    assert(pos >= 0 && pos < totLen);
    return fmasks[pos >> NSHIFT] & mask_codes[pos & MASK];
  }
  
  void setMask(int pos) {
    assert(pos >= 0 && pos < totLen);
    fmasks[pos >> NSHIFT] |= mask_codes[pos & MASK];
  }
  
private:
  int fullLen; // fullLen : the original length of an isoform
  int totLen; // totLen : the total length, included polyA tails, if any
  std::string name; // the tag
  std::string seq; // the raw sequence, in forward strand
  std::vector<uint32_t> fmasks; // record masks for forward strand, each position occupies 1 bit

  char dir;

  static const int NBITS = 32; // use unsigned int, 32 bits per variable
  static const int NSHIFT = 5;
  static const int MASK = (1 << NSHIFT) - 1;
  static const std::vector<uint32_t> mask_codes;

  static std::vector<uint32_t> init_mask_code();

  static const std::vector<int> base2code;
  static const std::vector<int> rbase2code;

  static std::vector<int> init_base2code();
  static std::vector<int> init_rbase2code();
};

#endif
