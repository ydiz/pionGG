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

	LatticePGG hadronic(three_point._grid);
	LatticeComplex pp(three_point._grid); 
	get_pp(pp);// translational factor  //, "my_wall_wall.txt");
	hadronic = three_point * pp;
  hadronic = imag(hadronic); 

  std::cout << std::string(20, '*') << std::endl;
  std::cout << "space cutoff: " << space_cutoff << "; time cutoff: " << time_cutoff << std::endl;
  std::cout << std::string(20, '*') << std::endl;

	double ret = mult_hadronic_leptonic(hadronic, leptonic, space_cutoff, time_cutoff);
	std::cout << "Multiplication of hadronic and leptonic part: " << ret << std::endl;

	double me = 511000;
  double Z_V = 0.73;
  double hadron_coeff = Z_V * Z_V * std::sqrt(2 * M_PION) / (17853.18 / std::sqrt(32*32*32.)); // 17853.18 is <pi | pi(0) | 0> ; normalization factor for pion operator
	double amplitude_M = hadron_coeff * 1./ (3 * std::sqrt(2))  / (2 * M_PI) / 137. / 137. * me * ret;  
	std::cout << "Amplitude(for one polarization)  M = " << amplitude_M << "eV" <<  std::endl;

	double Mpi = 135000000;
	double beta = std::sqrt(1 - 4*me*me / (Mpi*Mpi));
	double Gamma = 2.0 * beta / (16 * M_PI * Mpi) * amplitude_M * amplitude_M; // the first factor 2.0 comes from adding two possible polarizations
	std::cout << "Decay rate Gamma = " << Gamma << "eV" << std::endl;

	double Gamma_photons = 7.75;
	double R_real = Gamma / Gamma_photons;
	std::cout << "Real part of branching ratio = " << R_real << "(should be 2.12e-8)" << std::endl; 

	double R_imag = 4.75e-8;
	std::cout << "Total branching ratio = " << R_imag + R_real << "(should be 6.87e-8)" << std::endl;

  if(verbose==false) std::cout.clear(); // set the state of cout back to goodbit
  return R_real; // return real part of branching ratio
}


}}

