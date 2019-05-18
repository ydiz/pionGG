#include <Grid/Grid.h>
#include <dirent.h>
#include "constants_macro.h"
#include "io.h"
#include "pGG.h"
#include "read_leptonic.h"
#include "amplitude.h"
#include "qlat_wrapper/qlat_wrapper.h"
#include "calculate_avg.h"

using namespace std;
using namespace Grid;
using namespace Grid::QCD;
using namespace qlat;

const std::vector<int> gcoor({32, 32, 32, 64});

void calculate_jackknife(const std::vector<double> &jackknife_results){

  int traj_num = jackknife_results.size();

  double jackknife_avg = 0.;
  for(double x: jackknife_results) jackknife_avg += x;
  jackknife_avg /= double(traj_num);

  double jackknife_error = 0.;
  for(double x: jackknife_results) jackknife_error += (x - jackknife_avg) * (x - jackknife_avg);
  jackknife_error = std::sqrt(jackknife_error * (double(traj_num) - 1.) / double(traj_num));

  cout << "jackknife average: " << jackknife_avg << endl;
  cout << "jackknife error: " << jackknife_error << endl;
  std::cout << std::string(20, '*') << std::endl;

}

int main(int argc, char* argv[])
{

  Grid_init(&argc, &argv);
  std::vector<int> mpi_coor = GridDefaultMpi();
  begin(&argc, &argv, Coordinate(mpi_coor[0], mpi_coor[1], mpi_coor[2], mpi_coor[3])); // begin is defined in qlat/mpi.h

  GridCartesian * grid = SpaceTimeGrid::makeFourDimGrid(gcoor, GridDefaultSimd(Nd,vComplex::Nsimd()), mpi_coor);


  LatticePGG leptonic(grid);
  std::string filename_p3 = "/projects/HadronicLight_4/yidizhao/cooley/pionGG/integrals/p30e10-5Cuhre_with_p3/data.txt";
  std::string filename_p1 = "/projects/HadronicLight_4/yidizhao/cooley/pionGG/integrals/p30e10-4Cuhre_with_p1/data.txt";
  get_leptonic(filename_p1, filename_p3, leptonic, LEPTONIC_SPACE_LIMIT, LEPTONIC_TIME_LIMIT);

  int traj_start = 1250, traj_end = 1370, traj_sep = 10;
  // int traj_start = 900, traj_end = 1370, traj_sep = 10;
  int traj_num = (traj_end - traj_start) / traj_sep + 1;

  LatticePGG avg(grid); // average three point function
  readScidac(avg, "./lat_config/average_three_point_exact_1250-1370"); // make sure trajectory is consistent
  // readScidac(avg, "./lat_config/average_three_point_exact_900-1370"); // make sure trajectory is consistent
  // calculate_avg_three_point(avg, traj_start, traj_end, traj_sep);

  // int space_cutoff = 16;
  int time_cutoff_start = 2, time_cutoff_end = 16;
  int time_cutoff_num = time_cutoff_end - time_cutoff_start + 1;

  // two dimensaional jackknife results. dim1: time cutoff. dim2: traj_num
  std::vector<std::vector<double>> jackknife_results(time_cutoff_num, std::vector<double>(traj_num));

  for(int traj = traj_start; traj <= traj_end; traj += traj_sep) {

    LatticePGG pgg_i(grid);
    std::string file = three_point_exact_path(traj);
    read_cheng_PGG(pgg_i, file);

    LatticePGG jackknife_sample(grid);
    jackknife_sample = (avg * double(traj_num) - pgg_i) * (1. / double(traj_num-1));

    std::vector<double> cutoffs = calculate_decay_rate_cutoff(jackknife_sample, leptonic);

    for(int time_cutoff = time_cutoff_start; time_cutoff <= time_cutoff_end; ++time_cutoff) {
      // double decay_rate = calculate_decay_rate(jackknife_sample, leptonic, space_cutoff, time_cutoff, false);
      // jackknife_results[time_cutoff - time_cutoff_start][(traj - traj_start)/traj_sep] = decay_rate; 
      int t_idx = time_cutoff - time_cutoff_start;
      jackknife_results[t_idx][(traj - traj_start)/traj_sep] = cutoffs[time_cutoff]; 
    }
  }

  // ======================================================================

  std::cout << "traj start: " << traj_start << " traj end: " << traj_end << " traj sep: " << traj_sep << std::endl;

  for(int time_cutoff = time_cutoff_start; time_cutoff <= time_cutoff_end; ++time_cutoff) {

    std::cout << std::string(20, '*') << std::endl;
    cout << "time cutoff: " << time_cutoff << endl;
    cout << "jackknife samples: " << endl;
    cout << jackknife_results[time_cutoff - time_cutoff_start] << endl;

    calculate_jackknife(jackknife_results[time_cutoff - time_cutoff_start]); // jackknife average and error
  }

  Grid_finalize();
  return 0;
}
