#pragma once

namespace Grid {
namespace QCD {

void localIndexToLocalGlobalCoor(GridBase *grid, int ss, std::vector<int> &lcoor, std::vector<int> &gcoor) {
  // ss is local index; parallel_for(int ss=0; ss<ret._grid->lSites(); ss++)
  lcoor.resize(4);
  gcoor.resize(4);
  grid->LocalIndexToLocalCoor(ss, lcoor);
  std::vector<int> processor_coor;
  grid->ProcessorCoorFromRank(grid->ThisRank(), processor_coor);
  grid->ProcessorCoorLocalCoorToGlobalCoor(processor_coor, lcoor, gcoor);
}


// distance between two points with periodic boundary condition
int distance(int t1, int t2, int T) {
  int tmp = std::abs(t1 - t2);
  if(tmp <= T/2) return tmp;
  else return T-tmp;
}

int rightPoint(int t1, int t2, int T) {
  if(t1 > T/2 && t2 <= T/2) return t2;
  else if(t2 > T/2 && t1 <= T/2) return t1;
  else return std::max(t1, t2);
}


}}
