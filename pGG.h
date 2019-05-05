#pragma once

#include <qlat/qlat.h>
#include "io.h"
#include "constants_macro.h"

namespace Grid {
namespace QCD {

// to use Scidac, we make PGGElem has only one level (like iMatrix<Complex, 4>)
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

// get  exp(- 2 Mpi t) / sqrt( exp(Mpi t) <pi(0) pi(t)> ) for t in [0, T/2 - 10)
void get_pp(LatticeComplex &lat, const std::string &filename, int tmin=10) {

	int Tsize =  lat._grid->_fdimensions[Tdir];
	
	std::vector<double> pps = read_pp(filename);
	assert(pps.size() == Tsize);
	
	parallel_for(int ss=0; ss<lat._grid->lSites(); ss++){
		std::vector<int> lcoor(4);
		lat._grid->LocalIndexToLocalCoor(ss, lcoor);

		std::vector<int> gcoor(4);
		std::vector<int> processor_coor;
		lat._grid->ProcessorCoorFromRank(lat._grid->ThisRank(), processor_coor);
		lat._grid->ProcessorCoorLocalCoorToGlobalCoor(processor_coor, lcoor, gcoor);
		//std::cout << processor_coor << lcoor << gcoor << std::endl;

    // ===================================
    // cheng's tsep
    std::vector<int> x(4,0);
    std::vector<int> xp = gcoor;
    int t_min = 10;
    int t_wall;
    int t_sep;

    using namespace qlat;
    int diff = qlat::smod(x[3] - xp[3], 64);
    if (diff >= 0)
    {
      t_wall = qlat::mod(x[3] + t_min, 64);
    } else {
      t_wall = qlat::mod(xp[3] + t_min, 64);
    }
    t_sep = qlat::mod(t_wall - x[3], 64);

    // ===================================

		int xt = smod(gcoor[Tdir], Tsize);
		double val;
		if(xt <= TIME_LIMIT && xt >= -TIME_LIMIT) {
			int pion_t = (xt<=0) ? tmin : xt + tmin;
      assert(pion_t == t_wall); // FIXME: delete me
			val = 1. / std::sqrt(std::exp( - M_PION * pion_t) * pps[pion_t]); //FIXME: do I need to change this sign of pion_t?
			// val *= std::exp( - M_PION * pion_t); 
			// val *= std::exp( - M_PION * t_sep); // FIXME: cheng's t_sep is strange
		}
		else val = 0. ;
		
		typename LatticeComplex::vector_object::scalar_object m;
		m()()() = Complex(val, 0.);
		pokeLocalSite(m, lat, lcoor);
	}
}


}
}
