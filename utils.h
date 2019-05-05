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


}}
