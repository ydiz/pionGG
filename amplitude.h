#pragma once

#include <Grid/Grid.h>
#include <dirent.h>
#include "constants_macro.h"
#include "pGG.h"

namespace Grid{
namespace QCD{

std::vector<double> mult_HL_cutoff(const LatticePGG &hadronic, const LatticePGG &leptonic) {
  LatticeComplex tmp(hadronic._grid);
  tmp = 0.;

  parallel_for(int ss=0; ss<tmp._grid->oSites(); ++ss){
    tmp[ss]()()() = hadronic[ss]()()(0, 1) * leptonic[ss]()()(0, 1) + hadronic[ss]()()(1, 0) * leptonic[ss]()()(1, 0); 
    tmp[ss]()()() += hadronic[ss]()()(0, 2) * leptonic[ss]()()(0, 2) + hadronic[ss]()()(2, 0) * leptonic[ss]()()(2, 0); 
    tmp[ss]()()() += hadronic[ss]()()(2, 1) * leptonic[ss]()()(2, 1) + hadronic[ss]()()(1, 2) * leptonic[ss]()()(1, 2); 
  }

  // // int space_cutoff = 16;
  // int space_cutoff = 6;
  // std::cout << "space cutoff: " << space_cutoff << std::endl;
  // parallel_for(int ss=0; ss<tmp._grid->lSites(); ss++){
  //   std::vector<int> lcoor, gcoor;
  //   localIndexToLocalGlobalCoor(tmp._grid, ss, lcoor, gcoor);
  //
  //   if( (gcoor[0]>space_cutoff && gcoor[0]<tmp._grid->_fdimensions[0] - space_cutoff) ||
  //       (gcoor[1]>space_cutoff && gcoor[1]<tmp._grid->_fdimensions[1] - space_cutoff) ||
  //       (gcoor[2]>space_cutoff && gcoor[2]<tmp._grid->_fdimensions[2] - space_cutoff)  ) {
  //
  //     typename LatticeComplex::vector_object::scalar_object tmp_site; //  this has to be defined within parallel_for
  //     tmp_site = 0.;
  //     pokeLocalSite(tmp_site, tmp, lcoor);
  //   }
  // }

  int T = hadronic._grid->_fdimensions[Tdir];

  std::vector<iSinglet<Complex>> ret(T);
  sliceSum(tmp, ret, Tdir);

  std::vector<double> ret_real(T / 2 + 1);
  ret_real[0] = ret[0]()()().real();
  ret_real[T/2] = ret[T/2]()()().real();
  for(int i=1; i<ret_real.size()-1; ++i) ret_real[i] = ret[i]()()().real() + ret[T-i]()()().real(); // add t and -t

  std::vector<double> ret_cumulative(T / 2 + 1);
  ret_cumulative[0] = ret_real[0];
  for(int i=1; i<ret_real.size(); ++i) ret_cumulative[i] = ret_real[i] + ret_cumulative[i-1]; // result with cutoff
  return ret_cumulative;
}



// // L_{mu nu}(w) = lepton_coeff * leptonic
// // H_{mu nu}(w) = <0| Jmu(w/2) Jnu(-w/2) |pi> = hadron_coeff * three piont function
std::vector<double> calculate_decay_amplitude_cutoff(const LatticePGG &three_point, const LatticePGG &leptonic, double lepton_coeff, double hadron_coeff) { 
	std::vector<double> ret = mult_HL_cutoff(three_point, leptonic);

  std::vector<double> amplitude_M(ret.size());
  for(int i=0; i<ret.size(); ++i) amplitude_M[i] = hadron_coeff * lepton_coeff * ret[i];

  return amplitude_M; // return real part of branching ratio
}



}}

