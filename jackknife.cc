#include <Grid/Grid.h>
#include <dirent.h>
#include "constants_macro.h"
#include "io.h"
#include "pGG.h"
#include "read_leptonic.h"
#include "amplitude.h"
#include "qlat_wrapper/qlat_wrapper.h"
#include "jackknife.h"
#include "jack_init.h"

#include "imaginary_part.h"


using namespace std;
using namespace Grid;
using namespace Grid::QCD;
using namespace qlat;

const std::vector<int> gcoor({32, 32, 32, 64});

int main(int argc, char* argv[])
{

  Grid_init(&argc, &argv);
  std::vector<int> mpi_coor = GridDefaultMpi();
  begin(&argc, &argv, Coordinate(mpi_coor[0], mpi_coor[1], mpi_coor[2], mpi_coor[3])); // begin is defined in qlat/mpi.h

  Jack_para para;
  init_para(argc, argv, para);

  GridCartesian * grid = SpaceTimeGrid::makeFourDimGrid(gcoor, GridDefaultSimd(Nd,vComplex::Nsimd()), mpi_coor);

  // LatticePGG leptonic(grid);
  // std::string filename_p3 = "/projects/HadronicLight_4/yidizhao/cooley/pionGG/integrals/p30e10-5Cuhre_with_p3/data.txt";
  // std::string filename_p1 = "/projects/HadronicLight_4/yidizhao/cooley/pionGG/integrals/p30e10-4Cuhre_with_p1/data.txt";
  // get_leptonic(filename_p1, filename_p3, leptonic, LEPTONIC_SPACE_LIMIT, LEPTONIC_TIME_LIMIT);
  
  LatticeComplex leptonic(grid);
  imag_part(leptonic);
  //
  // cout << leptonic << endl;
  // assert(0);

  // two dimensaional jackknife results. dim1: time cutoff. dim2: traj_num
  std::vector<std::vector<double>> jackknife_results(para.time_cutoff_num, std::vector<double>(para.traj_num));

  for(int traj = para.traj_start; traj <= para.traj_end; traj += para.traj_sep) {

    LatticePGG three_point(grid);
    std::string file = three_point_exact_path(traj);
    read_cheng_PGG(three_point, file);

    // std::vector<double> cutoffs = calculate_decay_rate_cutoff(three_point, leptonic);

    std::vector<double> cutoffs = calculate_imag_decay_rate_cutoff(three_point, leptonic);

    for(int time_cutoff = para.time_cutoff_start; time_cutoff <= para.time_cutoff_end; ++time_cutoff) {
      int t_idx = time_cutoff - para.time_cutoff_start;
      jackknife_results[t_idx][(traj - para.traj_start)/para.traj_sep] = cutoffs[time_cutoff]; 
    }
  }

  // ======================================================================

  for(int time_cutoff = para.time_cutoff_start; time_cutoff <= para.time_cutoff_end; ++time_cutoff) {

    std::cout << std::string(20, '*') << std::endl;
    cout << "time cutoff: " << time_cutoff << endl;
    cout << "jackknife samples: " << endl;
    cout << jackknife_results[time_cutoff - para.time_cutoff_start] << endl;

    // calculate_jackknife(jackknife_results[time_cutoff - para.time_cutoff_start]); // jackknife average and error
    jack_stats(jackknife_results[time_cutoff - para.time_cutoff_start]); // jackknife average and error
  }

  Grid_finalize();
  return 0;
}
