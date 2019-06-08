#include <Grid/Grid.h>
#include <dirent.h>
#include "io.h"
#include "constants_macro.h"
#include "pGG.h"
#include "amplitude.h"
#include "qlat_wrapper/qlat_wrapper.h"
#include "jackknife.h"
#include "jack_init.h"
#include "lep.h"
#include "lep_CUBA3d.h"

#include "imaginary_part.h"


using namespace std;
using namespace Grid;
using namespace Grid::QCD;
using namespace qlat;

// const std::vector<int> gcoor({32, 32, 32, 64});
const std::vector<int> gcoor({24, 24, 24, 64});

int main(int argc, char* argv[])
{

  Grid_init(&argc, &argv);
  std::vector<int> mpi_coor = GridDefaultMpi();
  begin(&argc, &argv, Coordinate(mpi_coor[0], mpi_coor[1], mpi_coor[2], mpi_coor[3]));

  Jack_para para;
  init_para(argc, argv, para);

  GridCartesian * grid = SpaceTimeGrid::makeFourDimGrid(gcoor, GridDefaultSimd(Nd,vComplex::Nsimd()), mpi_coor);

  LatticePGG leptonic(grid);
  para.get_leptonic(leptonic); // generate leptonic part

  // two dimensaional jackknife results. dim1: time cutoff. dim2: traj_num
  std::vector<std::vector<double>> jackknife_results(para.time_cutoff_num, std::vector<double>(para.traj_num));

  for(int traj = para.traj_start; traj <= para.traj_end; traj += para.traj_sep) {

    LatticePGG three_point(grid);
    para.get_three_point(three_point, traj);

    std::vector<double> cutoffs = para.get_result_with_cutoff(three_point, leptonic);

    int traj_idx = (traj - para.traj_start)/para.traj_sep;
    for(int time_cutoff = para.time_cutoff_start; time_cutoff <= para.time_cutoff_end; ++time_cutoff) {
      int t_idx = time_cutoff - para.time_cutoff_start;
      jackknife_results[t_idx][traj_idx] = cutoffs[time_cutoff]; 
    }
  }

  // ======================================================================

  for(int time_cutoff = para.time_cutoff_start; time_cutoff <= para.time_cutoff_end; ++time_cutoff) {

    std::cout << std::string(20, '*') << std::endl;
    cout << "time cutoff: " << time_cutoff << endl;
    cout << "jackknife samples: " << endl;
    cout << jackknife_results[time_cutoff - para.time_cutoff_start] << endl;

    jack_stats(jackknife_results[time_cutoff - para.time_cutoff_start]); // jackknife average and error
  }

  Grid_finalize();
  return 0;
}
