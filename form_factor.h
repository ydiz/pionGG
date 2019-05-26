#pragma once

#include "utils.h"

namespace Grid {
namespace QCD {

void form_factor_integrand(LatticePGG &lat) {

	parallel_for(int ss=0; ss<lat._grid->lSites(); ss++){

    std::vector<int> lcoor, gcoor;
    localIndexToLocalGlobalCoor(lat._grid, ss, lcoor, gcoor);

    gcoor = my_smod(gcoor, lat._grid->_fdimensions);

    double w = std::sqrt(gcoor[0]*gcoor[0] + gcoor[1]*gcoor[1] + gcoor[2]*gcoor[2]);

    double val;
    if(w==0) val = 0.;
    else {
      double t = M_PION * 0.5 * w;
      double sin_t, cos_t;
      sincos(t, &sin_t, &cos_t);
      val = ( - cos_t * M_PION * w + 2 * sin_t) / w / w / w * std::exp(M_PION * 0.5 * gcoor[Tdir]);
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



std::vector<double> form_factor(const LatticePGG &three_point, const LatticePGG &leptonic) {
  // static LatticeComplex pp(three_point._grid); 
  // static bool pp_initialzed = false;
  // if(!pp_initialzed) {
  //   get_translational_factor(pp);// translational factor 
  //   pp_initialzed = true;
  // }

	LatticePGG hadronic(three_point._grid);
  hadronic = imag(three_point); 

	std::vector<double> ret = mult_HL_cutoff(hadronic, leptonic);

	double me = 511000;
  // double Z_V = 0.73;
  double Z_V = 0.7260;
  double hadron_coeff = 1./ (3 * std::sqrt(2)) * Z_V * Z_V * std::sqrt(2 * M_PION) / (17853.18 / std::sqrt(32*32*32.));
  double other_coeff = 8 * M_PI * M_PI * 0.092424 / M_PION / M_PION / M_PION / M_PION;
  std::vector<double> F(ret.size());
  for(int i=0; i<ret.size(); ++i) F[i] = hadron_coeff * other_coeff * ret[i];

  return F;
}



}}

