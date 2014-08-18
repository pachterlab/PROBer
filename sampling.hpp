#ifndef SAMPLING_H_
#define SAMPLING_H_

#include<ctime>
#include<cassert>
#include<vector>
#include<set>

#include<stdint.h>
#include "boost/random.hpp"

typedef uint32_t seedType;
typedef boost::random::mt19937 engine_type;
typedef boost::random::uniform_01<> uniform_01_dist;
typedef boost::random::gamma_distribution<> gamma_dist;
typedef boost::random::variate_generator<engine_type&, uniform_01_dist> uniform_01_generator;
typedef boost::random::variate_generator<engine_type&, gamma_dist> gamma_generator;

class Sampler {
public:
  Sampler(seedType seed) {
    engine = new engine_type(seed);
    rg = new uniform_01_generator(*engine, uniform_01_dist()); 
  }

  Sampler(engine_type *engine) {
    assert(engine != NULL);
    this->engine = engine;
    rg = new uniform_01_generator(*engine, uniform_01_dist()); 
  }

  ~Sampler() {
    delete engine;
    delete rg;
  }

  double random() { return (*rg)(); }  

  int sample(double*, int);
  int sample(std::vector<double>&, int);

private:
  engine_type *engine;
  uniform_01_generator *rg;
};

// arr should be cumulative!
// interval : [,)
// random number should be in [0, arr[len - 1])
// If by chance arr[len - 1] == 0.0, one possibility is to sample uniformly from 0...len-1
inline int Sampler::sample(double* arr, int len) {
  int l, r, mid;
  double prb = random() * arr[len - 1];
    
  l = 0; r = len - 1;
  while (l <= r) {
    mid = (l + r) / 2;
    if (arr[mid] <= prb) l = mid + 1;
    else r = mid - 1;
  }
  assert(l < len);
  
  return l;
}

inline int Sampler::sample(std::vector<double>& arr, int len) {
  int l, r, mid;
  double prb = random() * arr[len - 1];

  l = 0; r = len - 1;
  while (l <= r) {
    mid = (l + r) / 2;
    if (arr[mid] <= prb) l = mid + 1;
    else r = mid - 1;
  }
  assert(l < len);

  return l;
}

class EngineFactory {
public:
  EngineFactory() { seedEngine = NULL; }
  ~EngineFactory() { if (seedEngine != NULL) delete seedEngine; }

  void init(seedType seed = time(NULL)) { 
    seedEngine = new engine_type(seed); 
    seedSet.clear();
  }

  engine_type* new_engine() {
    seedType seed;

    do {
      seed = (*seedEngine)();
      iter = seedSet.find(seed);
    } while (iter != seedSet.end());
    seedSet.insert(seed);
    
    return new engine_type(seed);
  }

  Sampler* new_sampler() {
    return (new Sampler(new_engine()));
  }

private:
  engine_type *seedEngine;
  std::set<seedType> seedSet; // Empty set of seeds
  std::set<seedType>::iterator iter; 
};

#endif
