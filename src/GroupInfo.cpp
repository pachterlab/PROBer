#include<cstdio>
#include<string>
#include<vector>

#include "my_assert.h"
#include "GroupInfo.hpp"

void GroupInfo::load(const char* groupF) {
  FILE *fi = fopen(groupF, "r");
  int pos;
  
  general_assert(fi != NULL, "Cannot open " + cstrtos(groupF) + "! It may not exist.");
  
  starts.clear();
  while(fscanf(fi, "%d", &pos) == 1) {
    starts.push_back(pos);
  }
  fclose(fi);
  
  m = starts.size() - 1;
  gids = new int[starts.back()];
  for (int i = 0; i < m; i++) {
    for (int j = starts[i]; j < starts[i + 1]; j++) {
      gids[j] = i;
    }
  }
}
