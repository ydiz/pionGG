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
      ret = (w_vec[Ydir] / w) * 1. / t * (cos_t - sin_t / t);
    }

    if(gcoor[Ydir] == 0) ret = 0.; // when z=0, hadronic part should be 0


    typename LatticeComplex::vector_object::scalar_object lat_site; //  this has to be defined within parallel_for
    lat_site()()() = Complex(ret, 0.);

    pokeLocalSite(lat_site, lat, lcoor);
  }
}


}}

