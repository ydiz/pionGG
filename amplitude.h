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

  int T = hadronic._grid->_fdimensions[Tdir];

  std::vector<iSinglet<Complex>> ret(T);
  sliceSum(tmp, ret, Tdir);

  std::vector<double> ret_real(T / 2 + 1);
  ret_real[0] = ret[0]()()().real();
  for(int i=1; i<ret_real.size(); ++i) ret_real[i] = ret[i]()()().real() + ret[T-i]()()().real(); // add t and -t
  
  std::vector<double> ret_cumulative(T / 2 + 1);
  ret_cumulative[0] = ret_real[0];
  for(int i=1; i<ret_real.size(); ++i) ret_cumulative[i] = ret_real[i] + ret_cumulative[i-1]; // result with cutoff
  return ret_cumulative;
}



// L_{mu nu}(w) = lepton_coeff * leptonic
std::vector<double> calculate_decay_rate_cutoff(const LatticePGG &three_point, const LatticePGG &leptonic, double lepton_coeff) {

	static LatticeComplex pp(three_point._grid); 
  static bool pp_initialzed = false;
  if(!pp_initialzed) {
    get_translational_factor(pp);// translational factor 
    pp_initialzed = true;
  }

	LatticePGG hadronic(three_point._grid);
	hadronic = three_point * pp;
  hadronic = imag(hadronic); 

	std::vector<double> ret = mult_HL_cutoff(hadronic, leptonic);

	double me = 511000;
  // double Z_V = 0.73;
  double Z_V = 0.7260;
  double hadron_coeff = 1./ (3 * std::sqrt(2)) * Z_V * Z_V * std::sqrt(2 * M_PION) / (17853.18 / std::sqrt(32*32*32.)); // 17853.18 = sqrt(<pi(0) pi(t)> exp(Mpi * |t|)) = <pi | pi(0) | 0> / sqrt(2 Mpi); normalization factor for pion operator
  // double lepton_coeff = 1. / (2 * M_PI) / 137. / 137. * me;
  std::vector<double> amplitude_M(ret.size());
  for(int i=0; i<ret.size(); ++i) amplitude_M[i] = hadron_coeff * lepton_coeff * ret[i];

	double Mpi = 135000000;
	double beta = std::sqrt(1 - 4*me*me / (Mpi*Mpi));
  double Gamma_coeff = 2.0 * beta / (16 * M_PI * Mpi); // the first factor 2.0 comes from adding two possible polarizations
	// double Gamma = Gamma_coeff * amplitude_M * amplitude_M;

	// double Gamma_photons = 7.75;
	double Gamma_photons = 7.82;
	// double R_real = Gamma / Gamma_photons;

  std::vector<double> R_real(ret.size());
  for(int i=0; i<ret.size(); ++i) R_real[i] = (Gamma_coeff / Gamma_photons) * amplitude_M[i] *amplitude_M[i]; 
  
  return R_real; // return real part of branching ratio
}



}}

