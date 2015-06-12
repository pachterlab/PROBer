#ifndef UTILS_H_
#define UTILS_H_

#include<cmath>
#include<cctype>
#include<string>
#include<vector>
#include<stdint.h>

typedef uint64_t HIT_INT_TYPE;
typedef uint64_t READ_INT_TYPE;

const int STRLEN = 10005 ;
const double EPSILON = 1e-300;

const int MASK_LEN = 24; // the last MASK_LEN bp of a sequence cannot be aligned if poly(A) tail is added

const char channelStr[2][STRLEN] = {"minus", "plus"};

extern bool verbose; // show detail intermediate outputs

//inline bool isZero(double a) { return fabs(a) < 1e-8; }
//inline bool isLongZero(double a) { return fabs(a) < 1e-30; }

/*
  In our context, isXXZero is called only for non-negative values. Thus we omit the fabs function
 */
inline bool isZero(double a) { return a < 1e-8; }
inline bool isLongZero(double a) { return a < 1e-30; }

inline std::string cleanStr(const std::string& str) {
  int len = str.length();
  int fr, to;

  fr = 0;
  while (fr < len && isspace(str[fr])) ++fr;
  to = len - 1;
  while (to >= 0 && isspace(str[to])) --to;

  return (fr <= to ? str.substr(fr, to - fr + 1) : "");
}

static const char code2base[] = "ACGTN";

static std::vector<char> init_base2rbase() {
  std::vector<char> vec(128, -1);
  vec['a'] = 't'; vec['A'] = 'T';
  vec['c'] = 'g'; vec['C'] = 'G';
  vec['g'] = 'c'; vec['G'] = 'C';
  vec['t'] = 'a'; vec['T'] = 'A';
  vec['n'] = 'n'; vec['N'] = 'N';

  return vec;
}

static const std::vector<char> base2rbase = init_base2rbase();

static std::vector<int> init_base2code() {
  std::vector<int> vec(128, -1);
  vec['a'] = vec['A'] = 0;
  vec['c'] = vec['C'] = 1;
  vec['g'] = vec['G'] = 2;
  vec['t'] = vec['T'] = 3;
  vec['n'] = vec['N'] = 4;
  
  return vec;
}

static const std::vector<int> base2code = init_base2code();

static std::vector<int> init_rbase2code() {
  std::vector<int> vec(128, -1);
  vec['a'] = vec['A'] = 3;
  vec['c'] = vec['C'] = 2;
  vec['g'] = vec['G'] = 1;
  vec['t'] = vec['T'] = 0;
  vec['n'] = vec['N'] = 4;
  
  return vec;
}

static const std::vector<int> rbase2code = init_rbase2code();

#endif
