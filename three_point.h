#pragma once

#include <qlat/qlat.h>
#include <qlat/grid.h>
#include <dirent.h>
#include <sys/stat.h>
#include "qlat_wrapper/qlat_wrapper.h"
#include "constants_macro.h"
#include "io.h"
#include "utils.h"

namespace Grid {
namespace QCD{


void read_wall_src_props(const std::string &wall_src_path, const std::string &gauge_transform_path, std::vector<LatticePropagator> &wall_props) {
  using namespace qlat;
  assert(wall_props.size() == 64);

  // read gauge transformation
  std::cout << "Load Gauge Transform And Get Inv: " <<  gauge_transform_path << std::endl;
  assert(dirExists(gauge_transform_path));
  GaugeTransform qlat_gtinv;
  {
    GaugeTransform qlat_gt;
    dist_read_field(qlat_gt, gauge_transform_path);
    to_from_big_endian_64(get_data(qlat_gt)); 
    gt_inverse(qlat_gtinv, qlat_gt);
  }
  LatticeColourMatrix gt(wall_props[0]._grid);
  grid_convert(gt, qlat_gtinv);

  std::vector<int> ts;
  std::map<int, std::string> wall_subdirs;
  get_ts(wall_src_path, ts, wall_subdirs);
  assert(dirExists(wall_subdirs[0]));

  std::cout << "reading wall source propagators and applying gauge transformations" << std::endl;
  for(int t=0; t<wall_props[0]._grid->_fdimensions[Tdir]; ++t) {
    read_propagator(wall_props[t], wall_subdirs[t]);
    wall_props[t] = gt * wall_props[t];
  }
}


LatticePGG three_point_contraction(const std::vector<LatticePropagator> &wall_props, const LatticePropagator &point_prop, const std::vector<int> &xg) {
  LatticePGG ret(point_prop._grid);

  // !! global peekSite, not local peekLocalSite
  // !! cannot use global peekSite inside parallel_for; every node must peek the same site
  std::vector<typename LatticePropagator::vector_object::scalar_object> wall_props_to_xg(ret._grid->_fdimensions[Tdir]);
  for(int i=0; i<ret._grid->_fdimensions[Tdir]; ++i)	peekSite(wall_props_to_xg[i], wall_props[i], xg); 

	Gamma gamma5(Gamma::Algebra::Gamma5);
	// Gamma::gmu[0], Gamma::gmu[1], Gamma::gmu[0], Gamma::gmu[1];// GammaX, GammaY, GammaZ, GammaT

  parallel_for(int ss=0; ss<ret._grid->lSites(); ss++){

    std::vector<int> lcoor, gcoor;
    localIndexToLocalGlobalCoor(ret._grid, ss, lcoor, gcoor);

    std::vector<int> xp = gcoor;
    std::vector<int> x = xg;

    // cheng's tsep
    int t_min = 10;
    int t_wall;
    int t_sep;
    // using namespace qlat;
    int diff = qlat::smod(x[3] - xp[3], 64);
    if (diff >= 0)
    {
      t_wall = qlat::mod(x[3] + t_min, 64);
    } else {
      t_wall = qlat::mod(xp[3] + t_min, 64);
    }
    t_sep = qlat::mod(t_wall - x[3], 64);

    typename LatticePropagator::vector_object::scalar_object wall_to_x, x_to_wall, wall_to_xp, xp_to_wall, x_to_xp, xp_to_x;
    typename LatticePGG::vector_object::scalar_object ret_site;

    wall_to_x = wall_props_to_xg[t_wall];
    x_to_wall = gamma5 * adj(wall_to_x) * gamma5;
    peekLocalSite(wall_to_xp, wall_props[t_wall], lcoor);
    xp_to_wall = gamma5 * adj(wall_to_xp) * gamma5;
    peekLocalSite(x_to_xp, point_prop, lcoor);
    xp_to_x = gamma5 * adj(x_to_xp) * gamma5;

    ret_site = 0.;
    for(int mu=0; mu<4; ++mu)
      for(int nu=0; nu<4; ++nu) {
        ret_site()()(mu, nu) = trace(Gamma::gmu[mu] * xp_to_x * Gamma::gmu[nu] * wall_to_xp * gamma5 * x_to_wall); // cheng's order
        // ret_site()()(mu, nu) = trace(Gamma::gmu[mu] * xp_to_x * Gamma::gmu[nu] * wall_to_xp * gamma5 * x_to_wall) + trace(Gamma::gmu[nu] * x_to_xp * Gamma::gmu[mu] * wall_to_x * gamma5 * xp_to_wall); // cheng's order
    }

    ret_site = ret_site * (1 /  std::exp(- M_PION * t_sep));

    pokeLocalSite(ret_site, ret, lcoor);

  }
  return ret;
}

}}
