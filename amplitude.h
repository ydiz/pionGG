#pragma once

#include <Grid/Grid.h>
#include <dirent.h>
#include "constants_macro.h"
#include "io.h"
#include "pGG.h"
#include "read_leptonic.h"

namespace Grid{
namespace QCD{


bool isWithinRegion(const std::vector<int> gcoor, int space_cutoff, int time_cutoff, const std::vector<int> fdims) {
  if((gcoor[0] <= space_cutoff || gcoor[0] >= fdims[0] - space_cutoff) && (gcoor[1] <= space_cutoff || gcoor[1] >= fdims[1] - space_cutoff) && (gcoor[2] <= space_cutoff || gcoor[2] >= fdims[2] - space_cutoff) && (gcoor[3] <= time_cutoff || gcoor[3] >= fdims[3] - time_cutoff)) return true;
  else return false;
}


double mult_hadronic_leptonic(const LatticePGG &hadronic, const LatticePGG &leptonic, int space_cutoff, int time_cutoff) {
	LatticeComplex tmp(hadronic._grid);
  tmp = 0.;

  parallel_for(int ss=0; ss<tmp._grid->lSites(); ss++){
    std::vector<int> lcoor, gcoor;
    localIndexToLocalGlobalCoor(tmp._grid, ss, lcoor, gcoor);

    if(!isWithinRegion(gcoor, space_cutoff, time_cutoff, tmp._grid->_fdimensions)) continue; // set cutoff

    typename LatticeComplex::vector_object::scalar_object tmp_site; //  this has to be defined within parallel_for
		typename LatticePGG::vector_object::scalar_object hadronic_site, leptonic_site;
		peekLocalSite(hadronic_site, hadronic, lcoor);
		peekLocalSite(leptonic_site, leptonic, lcoor);

		tmp_site()()() = hadronic_site()()(0, 1) * leptonic_site()()(0, 1) + hadronic_site()()(1, 0) * leptonic_site()()(1, 0); 
		tmp_site()()() += hadronic_site()()(0, 2) * leptonic_site()()(0, 2) + hadronic_site()()(2, 0) * leptonic_site()()(2, 0); 
		tmp_site()()() += hadronic_site()()(2, 1) * leptonic_site()()(2, 1) + hadronic_site()()(1, 2) * leptonic_site()()(1, 2); 
    pokeLocalSite(tmp_site, tmp, lcoor);
  }

	Complex ret = TensorRemove(sum(tmp));
	return ret.real();
}


double calculate_decay_rate(const LatticePGG &three_point, const LatticePGG &leptonic, int space_cutoff, int time_cutoff, bool verbose=true) {

  if(verbose==false) std::cout.clear(std::ios_base::eofbit); // set state of cout to eofbit. cout will not print anything

	static LatticeComplex pp(three_point._grid); 
  static bool pp_initialzed = false;
  if(!pp_initialzed) {
    get_translational_factor(pp);// translational factor  //, "my_wall_wall.txt");
    pp_initialzed = true;
  }

	LatticePGG hadronic(three_point._grid);
	hadronic = three_point * pp;
  hadronic = imag(hadronic); 

  std::cout << std::string(20, '*') << std::endl;
  std::cout << "space cutoff: " << space_cutoff << "; time cutoff: " << time_cutoff << std::endl;
  std::cout << std::string(20, '*') << std::endl;

	double ret = mult_hadronic_leptonic(hadronic, leptonic, space_cutoff, time_cutoff);
	std::cout << "Multiplication of hadronic and leptonic part: " << ret << std::endl;

	double me = 511000;
  double Z_V = 0.73;
  double hadron_coeff = 1./ (3 * std::sqrt(2)) * Z_V * Z_V * std::sqrt(2 * M_PION) / (17853.18 / std::sqrt(32*32*32.)); // 17853.18 is <pi | pi(0) | 0> ; normalization factor for pion operator
  double lepton_coeff = 1. / (2 * M_PI) / 137. / 137. * me;
	double amplitude_M = hadron_coeff * lepton_coeff * ret;  
	std::cout << "Amplitude(for one polarization)  M = " << amplitude_M << "eV" <<  std::endl;

	double Mpi = 135000000;
	double beta = std::sqrt(1 - 4*me*me / (Mpi*Mpi));
  double Gamma_coeff = 2.0 * beta / (16 * M_PI * Mpi); // the first factor 2.0 comes from adding two possible polarizations
	double Gamma = Gamma_coeff * amplitude_M * amplitude_M;
	std::cout << "Decay rate Gamma = " << Gamma << "eV" << std::endl;

	double Gamma_photons = 7.75;
	double R_real = Gamma / Gamma_photons;
	std::cout << "Real part of branching ratio = " << R_real << "(should be 2.12e-8)" << std::endl; 

	double R_imag = 4.75e-8;
	std::cout << "Total branching ratio = " << R_imag + R_real << "(should be 6.87e-8)" << std::endl;

  if(verbose==false) std::cout.clear(); // set the state of cout back to goodbit
  return R_real; // return real part of branching ratio
}




// ============  for jackknife:

std::vector<double> mult_hadronic_leptonic_cutoff(const LatticePGG &hadronic, const LatticePGG &leptonic) {
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

std::vector<double> calculate_decay_rate_cutoff(const LatticePGG &three_point, const LatticePGG &leptonic) {

	static LatticeComplex pp(three_point._grid); 
  static bool pp_initialzed = false;
  if(!pp_initialzed) {
    get_translational_factor(pp);// translational factor 
    pp_initialzed = true;
  }

	LatticePGG hadronic(three_point._grid);
	hadronic = three_point * pp;
  hadronic = imag(hadronic); 

	std::vector<double> ret = mult_hadronic_leptonic_cutoff(hadronic, leptonic);

	double me = 511000;
  double Z_V = 0.73;
  double hadron_coeff = 1./ (3 * std::sqrt(2)) * Z_V * Z_V * std::sqrt(2 * M_PION) / (17853.18 / std::sqrt(32*32*32.)); // 17853.18 is <pi | pi(0) | 0> ; normalization factor for pion operator
  double lepton_coeff = 1. / (2 * M_PI) / 137. / 137. * me;
  std::vector<double> amplitude_M(ret.size());
  for(int i=0; i<ret.size(); ++i) amplitude_M[i] = hadron_coeff * lepton_coeff * ret[i];

	double Mpi = 135000000;
	double beta = std::sqrt(1 - 4*me*me / (Mpi*Mpi));
  double Gamma_coeff = 2.0 * beta / (16 * M_PI * Mpi); // the first factor 2.0 comes from adding two possible polarizations
	// double Gamma = Gamma_coeff * amplitude_M * amplitude_M;

	double Gamma_photons = 7.75;
	// double R_real = Gamma / Gamma_photons;

  std::vector<double> R_real(ret.size());
  for(int i=0; i<ret.size(); ++i) R_real[i] = (Gamma_coeff / Gamma_photons) * amplitude_M[i] *amplitude_M[i]; 
  
  return R_real; // return real part of branching ratio
}









// ===========================================================
// imaginary part

std::vector<double> imag_mult_hadronic_leptonic_cutoff(const LatticePGG &hadronic, const LatticeComplex &leptonic) {
	LatticeComplex tmp(hadronic._grid);
  // tmp = 0.;
  //
  // parallel_for(int ss=0; ss<tmp._grid->oSites(); ++ss){
	// 	tmp[ss]()()() = hadronic[ss]()()(0, 1) * leptonic[ss]()()() - hadronic[ss]()()(1, 0) * leptonic[ss]()()(); 
	// 	// tmp[ss]()()() += hadronic[ss]()()(0, 2) * leptonic[ss]()()(0, 2) + hadronic[ss]()()(2, 0) * leptonic[ss]()()(2, 0); 
	// 	// tmp[ss]()()() += hadronic[ss]()()(2, 1) * leptonic[ss]()()(2, 1) + hadronic[ss]()()(1, 2) * leptonic[ss]()()(1, 2); 
  // }
  
  // tmp = (peekColour(hadronic, Xdir, Ydir) - peekColour(hadronic, Ydir, Xdir)) * leptonic;
  // std::cout << tmp  << std::endl;
  // assert(0);
  tmp = (peekColour(hadronic, Xdir, Zdir) - peekColour(hadronic, Zdir, Xdir)) * leptonic;
  // tmp = (peekColour(hadronic, Ydir, Zdir) - peekColour(hadronic, Zdir, Ydir)) * leptonic;

  int T = hadronic._grid->_fdimensions[Tdir];

  std::vector<iSinglet<Complex>> ret(T);
  sliceSum(tmp, ret, Tdir);

  // std::cout << ret << std::endl;

  std::vector<double> ret_real(T / 2 + 1);
  ret_real[0] = ret[0]()()().real();
  for(int i=1; i<ret_real.size(); ++i) ret_real[i] = ret[i]()()().real() + ret[T-i]()()().real(); // add t and -t
  
  std::vector<double> ret_cumulative(T / 2 + 1);
  ret_cumulative[0] = ret_real[0];
  for(int i=1; i<ret_real.size(); ++i) ret_cumulative[i] = ret_real[i] + ret_cumulative[i-1]; // result with cutoff
  return ret_cumulative;
}




std::vector<double> calculate_imag_decay_rate_cutoff(const LatticePGG &three_point, const LatticeComplex &leptonic) {

	static LatticeComplex pp(three_point._grid); 
  static bool pp_initialzed = false;
  if(!pp_initialzed) {
    get_translational_factor(pp);// translational factor 
    pp_initialzed = true;
  }

	LatticePGG hadronic(three_point._grid);
	hadronic = three_point * pp;
  hadronic = imag(hadronic); 

    // print_grid_field_site(hadronic, std::vector<int>{1,2,3,4});
    // print_grid_field_site(leptonic, std::vector<int>{1,2,3,4});
	std::vector<double> ret = imag_mult_hadronic_leptonic_cutoff(hadronic, leptonic);

	double me = 511000;
	double Mpi = 135000000;
	double beta = std::sqrt(1 - 4*me*me / (Mpi*Mpi));
  double Z_V = 0.73;

  double hadron_coeff = 1./ (3 * std::sqrt(2)) * Z_V * Z_V * std::sqrt(2 * M_PION) / (17853.18 / std::sqrt(32*32*32.)); // 17853.18 is <pi | pi(0) | 0> ; normalization factor for pion operator // FIXME: space volume and phi ground state amplitude
  // double lepton_coeff = me / Mpi * M_PI / 137. / 137. * (1. / beta * std::log((1 + beta) / (1 - beta)));
  double lepton_coeff = me / M_PION * M_PI / 137. / 137. * (1. / beta * std::log((1 + beta) / (1 - beta)));

  std::vector<double> amplitude_M(ret.size());
  // for(int i=0; i<ret.size(); ++i) amplitude_M[i] = hadron_coeff * lepton_coeff * ret[i];
  for(int i=0; i<ret.size(); ++i) amplitude_M[i] = 3.0 * hadron_coeff * lepton_coeff * ret[i]; //FIXME: 3.0

  double Gamma_coeff = 2.0 * beta / (16 * M_PI * Mpi); // the first factor 2.0 comes from adding two possible polarizations
	// double Gamma = Gamma_coeff * amplitude_M * amplitude_M;

	double Gamma_photons = 7.75;
	// double R_real = Gamma / Gamma_photons;

  std::vector<double> R_real(ret.size());
  for(int i=0; i<ret.size(); ++i) R_real[i] = (Gamma_coeff / Gamma_photons) * amplitude_M[i] *amplitude_M[i]; 
  
  // return R_real; // return real part of branching ratio
  return R_real; // return real part of branching ratio
}



}}

