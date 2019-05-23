#pragma once

#include "utils.h"


namespace Grid {
namespace QCD {



void imag_part(LatticeComplex &lat) {
  std::vector<int> space_vol(lat._grid->_fdimensions.begin(), lat._grid->_fdimensions.begin()+3);

  parallel_for(int ss=0; ss<lat._grid->lSites(); ss++){
    std::vector<int> lcoor, gcoor;
    localIndexToLocalGlobalCoor(lat._grid, ss, lcoor, gcoor);

    std::vector<int> space_coor(gcoor.begin(), gcoor.begin()+3);

    double ret;

    std::vector<int> w_vec = my_smod(space_coor, space_vol);

    double w = len(w_vec);

    if(w==0.) ret = 0.;
    else {
      double t = M_PION * 0.5 * w;
      double sin_t, cos_t;
      sincos(t, &sin_t, &cos_t);
      ret = (w_vec[Zdir] / w) * 1. / t * (cos_t - sin_t / t);
    }

    // if(gcoor[3]>16 || gcoor[3] <64-16) ret = 0;


    typename LatticeComplex::vector_object::scalar_object lat_site; //  this has to be defined within parallel_for
    lat_site()()() = Complex(ret, 0.);

    pokeLocalSite(lat_site, lat, lcoor);
  }
}


void imag_part(LatticePGG &lat) {

	parallel_for(int ss=0; ss<lat._grid->lSites(); ss++){

    std::vector<int> lcoor, gcoor;
    localIndexToLocalGlobalCoor(lat._grid, ss, lcoor, gcoor);

    gcoor = my_smod(gcoor, lat._grid->_fdimensions);

    double w = std::sqrt(gcoor[0]*gcoor[0] + gcoor[1]*gcoor[1] + gcoor[2]*gcoor[2] );

		// double val = get_integral_site(gcoor, data, lat._grid->_fdimensions, space_limit, time_limit); 
    double val;
    if(w==0) val = 0.;
    else {
      double t = M_PION * 0.5 * w;
      double sin_t, cos_t;
      sincos(t, &sin_t, &cos_t);
      val = (1. / w) * 1. / t * (cos_t - sin_t / t);
    }

		typename LatticePGG::vector_object::scalar_object m;
		m = 0.;
		m()()(0, 1) = Complex(val * gcoor[Zdir], 0); 
		m()()(0, 2) = Complex(-val * gcoor[Ydir], 0); // Minus sign comes from spinor matrix
		m()()(1, 2) = Complex(val * gcoor[Xdir], 0); 
		m()()(1, 0) = - m()()(0, 1);
		m()()(2, 0) = - m()()(0, 2);
		m()()(2, 1) = - m()()(1, 2);

		pokeLocalSite(m, lat, lcoor);
	}
}


}}

