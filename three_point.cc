#include "three_point.h"

std::vector<int> gcoor({32, 32, 32, 64});
std::vector<int> mpi_coor({1, 1, 1, 8});

using namespace std;
using namespace qlat;
using namespace Grid;
using namespace Grid::QCD;


int main(int argc, char* argv[])
{
  begin(&argc, &argv, Coordinate(mpi_coor[0], mpi_coor[1], mpi_coor[2], mpi_coor[3])); // begin is defined in qlat/mpi.h
  Grid_init(&argc, &argv);

  GridCartesian *grid = SpaceTimeGrid::makeFourDimGrid(gcoor, GridDefaultSimd(Nd,vComplex::Nsimd()), mpi_coor); 

	// int traj_start = 680, traj_num = 70;
	// std::vector<int> trajs(traj_num);
	// for(int i=0; i<trajs.size(); ++i) trajs[i] = traj_start + i * 10;
	std::vector<int> trajs {1370};

	cout << "trajs: " << endl;
	cout << trajs << endl;

	for(int traj: trajs) {

		std::string gauge_transform_path = gauge_transform_path_32D(traj);
		std::string wall_src_path = wall_path_32D_sloppy(traj);
		std::string point_src_path = point_path_32D(traj);

    // read wall src
		std::vector<LatticePropagator> wall_props(gcoor[Tdir], grid);
    read_wall_src_props(wall_src_path, gauge_transform_path, wall_props);

		std::vector<std::vector<int>> xgs;
		std::map<std::vector<int>, std::string> point_subdirs;
		get_xgs(point_src_path, xgs, point_subdirs);
		xgs.clear(); xgs.push_back({0, 4, 2, 6});

		for(const auto &xg: xgs) {

			cout << "xg of point src: " << xg << endl;
			LatticePropagator point_prop(grid);
			read_propagator(point_prop, point_subdirs[xg]);
			std::cout << "Finished reading point propagator!" << std::endl;

      LatticePGG ret = three_point_contraction(wall_props, point_prop, xg);

			for(int mu=0; mu<4; ++mu) ret = Cshift(ret, mu, xg[mu]);

      ret = real(ret);
      ret = 2.0 * ret; // two diagrams: clockwise and anti-clockwise; they are conjugate to each other.

      // writeScidac(ret, "xg=(0,4,2,6)");
			cout << ret << endl;
		}

	}


  end();

  return 0;
}
