#pragma once

#include <Grid/Grid.h>
#include <qlat/qlat.h>
#include <dirent.h>
#include "io.h"
#include "pGG.h"
#include "qlat_wrapper/qlat_wrapper.h"
#include "constants_macro.h"


namespace Grid {
namespace QCD {

void calculate_avg_three_point(LatticePGG &avg, int traj_start, int traj_end, int traj_sep) {

  avg = 0.;

  int traj_num = (traj_end - traj_start) / traj_sep + 1;
	std::vector<int> trajs(traj_num);
	for(int i=0; i<trajs.size(); ++i) trajs[i] = traj_start + i * traj_sep;

  std::cout << "trajs: " << std::endl;
	std::cout << trajs << std::endl;

	for(int traj: trajs) {

    LatticePGG pgg(avg._grid);
    // std::string file = three_point_sloppy_path(traj);
    std::string file = three_point_exact_path(traj);
    read_cheng_PGG(pgg, file);

    avg += pgg;
    // print_grid_field_site(pgg, std::vector<int>{1, 2, 3, 4});
	}

	avg = avg * (1. /  double(trajs.size()));
}




}}
