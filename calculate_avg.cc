#include <Grid/Grid.h>
// #include <qlat/qlat.h>
// #include <qlat/grid.h>
#include <dirent.h>
#include "io.h"
#include "pGG.h"
#include "qlat_wrapper/qlat_wrapper.h"
#include "constants_macro.h"


using namespace std;
using namespace qlat;
using namespace Grid;
using namespace Grid::QCD;

// void get_trajs(const std::string &path, std::vector<int> &trajs, std::map<int, std::string> &files) {
//     DIR *dir;
//     dir = opendir(path.c_str());
//     assert(dir!=NULL); // make sure directory exists
//     struct dirent *entry;
//
//     std::string filename;
//     while ((entry = readdir (dir)) != NULL) {
// 			// printf ("%s\n", entry->d_name);
// 			filename = std::string(entry->d_name);
// 			if(filename.substr(0, 4) == "PGG.") {
// 				int traj = std::stoi(filename.substr( filename.find(".") + 1 )); // "PGG.1000";
// 				trajs.push_back(traj);
// 				files.insert(std::pair<int, std::string>(traj, path + "/" + filename));
// 			}
//     }
//     closedir (dir);
//
//     std::sort(trajs.begin(), trajs.end());
// }


const std::vector<int> gcoor({32, 32, 32, 64});
const std::vector<int> mpi_coor {1,1,1,8};

int main(int argc, char* argv[])
{

  begin(&argc, &argv, Coordinate(mpi_coor[0], mpi_coor[1], mpi_coor[2], mpi_coor[3])); // begin is defined in qlat/mpi.h
  Grid_init(&argc, &argv);

  GridCartesian * grid = SpaceTimeGrid::makeFourDimGrid(gcoor, GridDefaultSimd(Nd,vComplex::Nsimd()), mpi_coor);
	LatticePGG pgg(grid);
	LatticePGG avg(grid);
	avg = 0.;

	// std::string dir = "/home/ydzhao/cuth/three_point_config/32D_ama";

	// int traj_start = 690, traj_num = 69;
	int traj_start = 830, traj_num = 55;
	std::vector<int> trajs(traj_num);
	for(int i=0; i<trajs.size(); ++i) trajs[i] = traj_start + i * 10;

	cout << "trajs: " << endl;
	cout << trajs << endl;

	// std::vector<int> trajs; 
	// std::map<int, std::string> files;
	// get_trajs(dir, trajs, files);

	for(int traj: trajs) {

		// std::string file = three_point_sloppy_path(traj);
		std::string file = three_point_exact_path(traj);
		if(!dirExists(file)) cout << "The following file does not exist and I am skipping it: " << file << endl; 

		qlat::PionGGElemField qlat_pgg;

		cout << "reading from: " << file << endl;

		// readScidac(pgg, file);
		dist_read_field(qlat_pgg, file);
		cout << "Finished reading" << endl;
		grid_convert(pgg, qlat_pgg);

		avg += pgg;
		// print_grid_field_site(pgg, std::vector<int>{1, 13, 23, 47});
		// print_grid_field_site(pgg, std::vector<int>{1, 13, 23, 10});
		print_grid_field_site(pgg, std::vector<int>{0,0,0,0});
		print_grid_field_site(pgg, std::vector<int>{1, 2, 3, 4});
	}

	avg = avg * (1. /  double(trajs.size()));
	// print_grid_field_site(avg, std::vector<int>{1, 13, 23, 10});
	print_grid_field_site(avg, std::vector<int>{0,0,0,0});
	print_grid_field_site(avg, std::vector<int>{1, 2, 3, 4});

	writeScidac(avg, "./average_three_point_exact");
	
	// readScidac(avg, dir + "average");
	// LatticeComplex pp(grid); 
	// get_pp(pp, "wall_wall.txt"); // pp = 1 / sqrt( exp(Mpi t) <pi(0) pi(t)> ) for t in [0, T/2 - 10)
  // //
	// // matrix element = e^(-Mpi t) * three point / sqrt( exp(Mpi t) <pi(0) pi(t)> )
	// LatticePGG matrix_element(grid);
	// matrix_element = avg * pp;
	// writeScidac(matrix_element, dir + "/matrix_element");

	Grid_finalize();
  return 0;
}
