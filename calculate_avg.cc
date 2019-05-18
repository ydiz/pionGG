#include <Grid/Grid.h>
#include <qlat/qlat.h>
#include <dirent.h>
#include "io.h"
#include "pGG.h"
#include "qlat_wrapper/qlat_wrapper.h"
#include "constants_macro.h"
#include "calculate_avg.h"

using namespace std;
using namespace qlat;
using namespace Grid;
using namespace Grid::QCD;

const std::vector<int> gcoor({32, 32, 32, 64});

int main(int argc, char* argv[])
{
  Grid_init(&argc, &argv);
  std::vector<int> mpi_coor = GridDefaultMpi();
  begin(&argc, &argv, Coordinate(mpi_coor[0], mpi_coor[1], mpi_coor[2], mpi_coor[3])); // begin is defined in qlat/mpi.h

  GridCartesian * grid = SpaceTimeGrid::makeFourDimGrid(gcoor, GridDefaultSimd(Nd,vComplex::Nsimd()), mpi_coor);

	// int traj_start = 1250, traj_end = 1370, traj_sep = 10;
	int traj_start = 900, traj_end = 1370, traj_sep = 10;

  LatticePGG avg(grid);
  calculate_avg_three_point(avg, traj_start, traj_end, traj_sep);

	writeScidac(avg, "./average_three_point_exact_" + std::to_string(traj_start) + "-" + std::to_string(traj_end));

	Grid_finalize();
  return 0;
}
