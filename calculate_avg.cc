#include <Grid/Grid.h>
#include <qlat/qlat.h>
#include <dirent.h>
#include "io.h"
#include "pGG.h"
#include "qlat_wrapper/qlat_wrapper.h"
#include "constants_macro.h"


using namespace std;
using namespace qlat;
using namespace Grid;
using namespace Grid::QCD;


const std::vector<int> gcoor({32, 32, 32, 64});
// const std::vector<int> mpi_coor {1,1,1,8};

int main(int argc, char* argv[])
{
  Grid_init(&argc, &argv);
  std::vector<int> mpi_coor = GridDefaultMpi();
  begin(&argc, &argv, Coordinate(mpi_coor[0], mpi_coor[1], mpi_coor[2], mpi_coor[3])); // begin is defined in qlat/mpi.h

  GridCartesian * grid = SpaceTimeGrid::makeFourDimGrid(gcoor, GridDefaultSimd(Nd,vComplex::Nsimd()), mpi_coor);
	LatticePGG pgg(grid);
	LatticePGG avg(grid);
	avg = 0.;

	// // std::string dir = "/home/ydzhao/cuth/three_point_config/32D_ama";

	int traj_start = 1250, traj_end = 1370, traj_sep = 10;
  int traj_num = (traj_end - traj_start) / traj_sep + 1;
	std::vector<int> trajs(traj_num);
	for(int i=0; i<trajs.size(); ++i) trajs[i] = traj_start + i * traj_sep;

	cout << "trajs: " << endl;
	cout << trajs << endl;

	for(int traj: trajs) {

		// std::string file = three_point_sloppy_path(traj);
		std::string file = three_point_exact_path(traj);
		if(!dirExists(file)) cout << "The following file does not exist and I am skipping it: " << file << endl; 

		qlat::PionGGElemField qlat_pgg;

		cout << "reading from: " << file << endl;

		dist_read_field(qlat_pgg, file);
		cout << "Finished reading" << endl;
		grid_convert(pgg, qlat_pgg);

		avg += pgg;

		print_grid_field_site(pgg, std::vector<int>{0,0,0,0});
		print_grid_field_site(pgg, std::vector<int>{1, 2, 3, 4});
	}

	avg = avg * (1. /  double(trajs.size()));
	print_grid_field_site(avg, std::vector<int>{0,0,0,0});
	print_grid_field_site(avg, std::vector<int>{1, 2, 3, 4});

	// writeScidac(avg, "./average_three_point_exact");

	Grid_finalize();
  return 0;
}
