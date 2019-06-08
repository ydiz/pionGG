#include "../qlat_wrapper/qlat_wrapper.h"
// #include "three_point.h"

std::vector<int> get_xg(const std::string &path) {

	std::stringstream ss;
	ss.str(path.substr(path.find("(") + 1));

	std::vector<int> ret(4);
	for(int &x: ret) { 
		ss >> x; 
		ss.ignore(); // extract comma and ignore it
	}
	return ret;
}


void get_xgs(const std::string &path, std::vector<std::vector<int>> &xgs, std::map<std::vector<int>, std::string> &subdirs) {
	DIR *dir;
	dir = opendir(path.c_str());
	assert(dir!=NULL); // make sure directory exists
	struct dirent *entry;

	std::string subdir_name;
	while ((entry = readdir (dir)) != NULL) {
		// printf ("%s\n", entry->d_name);
		subdir_name = std::string(entry->d_name);
		if(subdir_name.substr(0, 3) == "xg=" && subdir_name.substr(subdir_name.find("type"), 6) == "type=0" && subdir_name.substr(subdir_name.find("accuracy"), 10) == "accuracy=0") {
			std::vector<int> xg = get_xg(subdir_name); 
			xgs.push_back(xg);
			subdirs.insert(std::pair<std::vector<int>, std::string>(xg, path + "/" + subdir_name));
		}
	}
	closedir (dir);
}




std::vector<int> gcoor({32, 32, 32, 64});

using namespace std;
using namespace qlat;
using namespace Grid;
using namespace Grid::QCD;


int main(int argc, char* argv[])
{
  Grid_init(&argc, &argv);
  std::vector<int> mpi_coor = GridDefaultMpi();
  begin(&argc, &argv, Coordinate(mpi_coor[0], mpi_coor[1], mpi_coor[2], mpi_coor[3]));
  GridCartesian * grid = SpaceTimeGrid::makeFourDimGrid(gcoor, GridDefaultSimd(Nd,vComplex::Nsimd()), mpi_coor);


	// int traj_start = 680, traj_num = 70;
	// std::vector<int> trajs(traj_num);
	// for(int i=0; i<trajs.size(); ++i) trajs[i] = traj_start + i * 10;
	std::vector<int> trajs {1370};

	cout << "trajs: " << endl;
	cout << trajs << endl;

	for(int traj: trajs) {
    LatticePropagator average(grid);

		// std::string gauge_transform_path = gauge_transform_path_32D(traj);
		// std::string wall_src_path = wall_path_32D_sloppy(traj);
    //
    // // read wall src
		// std::vector<LatticePropagator> wall_props(gcoor[Tdir], grid);
    // read_wall_src_props(wall_src_path, gauge_transform_path, wall_props);

		std::vector<std::vector<int>> xgs;
		std::map<std::vector<int>, std::string> point_subdirs;
		std::string point_src_path = point_path_32D(traj);
		get_xgs(point_src_path, xgs, point_subdirs);
		// xgs.clear(); xgs.push_back({0, 4, 2, 6});

		for(const auto &xg: xgs) {

			cout << "source point: " << xg << endl;
			cout << "path: " << point_subdirs[xg] << endl;
			// LatticePropagator point_prop(grid);
			// read_qlat_propagator(point_prop, point_subdirs[xg]);
			// std::cout << "Finished reading point propagator!" << std::endl;
      //
			// for(int mu=0; mu<4; ++mu) ret = Cshift(ret, mu, xg[mu]);
		}

	}


  end();

  return 0;
}
