#pragma once

#include <Grid/Grid.h>
#include "pGG.h"
#include "lep.h"
#include "lep_CUBA3d.h"
#include "imaginary_part.h"

namespace Grid{
namespace QCD{

struct Lep_para {
  std::string target; // real, real_CUBA3d, imag_analytic, imag_CUBA3d
  std::string file_p3;
  std::string file_p1;
  double lep_coef();
  void get_leptonic(LatticePGG &lat);
};

double Lep_para::lep_coef() {

  double me = 511000; // Unit is eV
  double ret;

  if(target=="real") ret = 2. / (M_PI) / 137. / 137. * me;
  else if(target=="real_CUBA3d") ret = 1. / (2 * M_PI) / 137. / 137. * me;
  else if(target=="imag_analytic") {
    double Mpi = 135000000;
    double beta = std::sqrt(1 - 4*me*me / (Mpi*Mpi));
    ret = me / M_PION * M_PI / 137. / 137. * (1. / beta * std::log((1 + beta) / (1 - beta)));
  }
  else if(target=="imag_CUBA3d") ret = 1. / 2. / 137. / 137. * me;
  else assert(0);

  return ret;
}


void Lep_para::get_leptonic(LatticePGG &lat) {
  if(target == "real_CUBA3d" || target == "imag_CUBA3d") get_leptonic_CUBA3d(file_p1, file_p3, lat, LEPTONIC_SPACE_LIMIT, LEPTONIC_TIME_LIMIT);
  else if(target == "real") Grid::QCD::get_leptonic(file_p3, lat, LEPTONIC_SPACE_LIMIT, LEPTONIC_TIME_LIMIT);
  else if(target == "imag_analytic") imag_part(lat);
  else assert(0);
}


}}
