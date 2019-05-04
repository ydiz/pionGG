#include <Grid/Grid.h>
#include <dirent.h>
#include "constants_macro.h"
#include "io.h"
#include "pGG.h"

using namespace std;
using namespace Grid;
using namespace Grid::QCD;

namespace Grid{
namespace QCD{

double mult_hadronic_leptonic(const LatticePGG &hadronic, const LatticePGG &leptonic) {
	LatticeComplex tmp(hadronic._grid);
	parallel_for(int ss=0; ss<hadronic._grid->oSites(); ++ss){
		// tmp[ss]()()() = 2. * ( hadronic[ss]()()(0, 1) * leptonic[ss]()()(0, 1) + hadronic[ss]()()(2, 2) * leptonic[ss]()()(2, 3) ) ; 
		tmp[ss]()()() = 2. * hadronic[ss]()()(0, 1) * leptonic[ss]()()(0, 1); 
	}

  tmp = abs(tmp);
  // std::cout<< tmp << std::endl;

	Complex ret = TensorRemove(sum(tmp));
	return ret.real();
}


double mult_hadronic_leptonic_print(const LatticePGG &hadronic, const LatticePGG &leptonic) {
	LatticeComplex tmp(hadronic._grid);
	parallel_for(int ss=0; ss<hadronic._grid->oSites(); ++ss){
		// tmp[ss]()()() = 2. * ( hadronic[ss]()()(0, 1) * leptonic[ss]()()(0, 1) + hadronic[ss]()()(2, 2) * leptonic[ss]()()(2, 3) ) ; 
		tmp[ss]()()() = 2. * hadronic[ss]()()(0, 1) * leptonic[ss]()()(0, 1); 
	}

  std::cout<< tmp << std::endl;

	Complex ret = TensorRemove(sum(tmp));
	return ret.real();
}

}}

const std::vector<int> gcoor({32, 32, 32, 64});

int main(int argc, char* argv[])
{

  Grid_init(&argc, &argv);

  GridCartesian * grid = SpaceTimeGrid::makeFourDimGrid(gcoor, GridDefaultSimd(Nd,vComplex::Nsimd()), GridDefaultMpi());
	LatticePGG hadronic(grid);
	LatticePGG leptonic(grid);

	// std::string hadronic_file = "/home/ydzhao/cuth/three_point_config/32D_ama/matrix_element";
	// std::string leptonic_file = "/home/ydzhao/cuth/pionGG/lat_leptonic";
	std::string hadronic_file = "./lat_config/matrix_element";
	std::string leptonic_file = "./lat_config/lat_leptonic";

	readScidac(hadronic, hadronic_file);
	// hadronic = timesI(hadronic); //this is done before writing matrix_elment // cheng's three piont function is imaginary
	// cout << hadronic << endl;
  //
	readScidac(leptonic, leptonic_file);

	double ret = mult_hadronic_leptonic(hadronic, leptonic);
	cout << "Multiplication of hadronic and leptonic part: " << ret << endl;

	double me = 511000;
	double amplitude_M = 1./ (3 * std::sqrt(2))  / (2 * M_PI) / 137 / 137 * me * ret; 
	cout << "Amplitude M = " << amplitude_M << "eV" <<  endl;

	double Mpi = 135000000;
	double beta = std::sqrt(1 - 4*me*me / (Mpi*Mpi));
	double Gamma = beta / (16 * M_PI * Mpi) * amplitude_M * amplitude_M;
	cout << "Decay rate Gamma = " << Gamma << "eV" << endl;

	double Gamma_photons = 7.75;
	double R_real = Gamma / Gamma_photons;
	cout << "Real part of branching ratio = " << R_real << "(should be 2.12e-8 or 2.82e-8)" << endl; 

	double R_imag = 4.75e-8;
	cout << "Total branching ratio = " << R_imag + R_real << "(should be 6.87e-8 or 7.57e-8)" << endl;

	mult_hadronic_leptonic_print(hadronic, leptonic);

	Grid_finalize();
  return 0;
}
