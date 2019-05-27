#pragma once

#include <qlat/qlat.h>
#include "io.h"
#include "constants_macro.h"
#include "utils.h"


namespace Grid {
namespace QCD {

// to use Scidac, we cannot make PGGElem have only one level (like iMatrix<Complex, 4>)
using PGGElem = iScalar<iScalar<iMatrix<Complex, 4>>>;
using vPGGElem = iScalar<iScalar<iMatrix<vComplex, 4>>>;
using LatticePGG = Lattice<vPGGElem>;

template<class T>
void print_grid_field_site(const T &field, const std::vector<int> coor) {
	using namespace Grid;
	std::cout << "[ " << coor[0] << " " << coor[1] << " " << coor[2] << " " << coor[3] << " ]" << std::endl;
	typename T::vector_object::scalar_object site;
	peekSite(site, field, coor);
	std::cout << site << std::endl;
}

// translational factor
void get_translational_factor(LatticeComplex &lat, double Mpi) {

	parallel_for(int ss=0; ss<lat._grid->lSites(); ss++){
    std::vector<int> lcoor, gcoor;
    localIndexToLocalGlobalCoor(lat._grid, ss, lcoor, gcoor);

		double val;
		int xt = qlat::smod(gcoor[Tdir], lat._grid->_fdimensions[Tdir]);
    val = std::exp( 0.5 * Mpi * xt); // translation factor

		typename LatticeComplex::vector_object::scalar_object m;
		m()()() = Complex(val, 0.);
		pokeLocalSite(m, lat, lcoor);
	}
}


}
}
