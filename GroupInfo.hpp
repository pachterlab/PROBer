#ifndef GROUPINFO_H_
#define GROUPINFO_H_

#include<cassert>
#include<vector>

class GroupInfo {
public:
  GroupInfo() { m = 0; starts.clear(); gids = NULL; }
  ~GroupInfo() { m = 0; starts.clear(); if (gids != NULL) delete[] gids; }
  
  void load(const char*);
  
  int getm() const { return m; }
  
  int gidAt(int sid) const {
    assert(sid > 0 && sid < starts.back());
    return gids[sid];
  }
  
  // sp : start position
  int spAt(int gid) const {
    assert(gid >= 0 && gid <= m);
    return starts[gid];
  }
  
private:
  int m; // m genes
  std::vector<int> starts; // genes' start positions
  int *gids; // hash
};

#endif /* GROUPINFO_H_ */
