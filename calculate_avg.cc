#include <Grid/Grid.h>
#include <dirent.h>
#include "constants_macro.h"
#include "io.h"
#include "pGG.h"

using namespace std;
using namespace Grid;
using namespace Grid::QCD;

void get_trajs(const std::string &path, std::vector<int> &trajs, std::map<int, std::string> &files) {
    DIR *dir;
    dir = opendir(path.c_str());
    assert(dir!=NULL); // make sure directory exists
    struct dirent *entry;

    std::string filename;
    while ((entry = readdir (dir)) != NULL) {
			// printf ("%s\n", entry->d_name);
			filename = std::string(entry->d_name);
			if(filename.substr(0, 4) == "PGG.") {
				int traj = std::stoi(filename.substr( filename.find(".") + 1 )); // "PGG.1000";
				trajs.push_back(traj);
				files.insert(std::pair<int, std::string>(traj, path + "/" + filename));
			}
    }
    closedir (dir);

    std::sort(trajs.begin(), trajs.end());
}


const std::vector<int> gcoor({32, 32, 32, 64});

int main(int argc, char* argv[])
{

  Grid_init(&argc, &argv);

  GridCartesian * grid = SpaceTimeGrid::makeFourDimGrid(gcoor, GridDefaultSimd(Nd,vComplex::Nsimd()), GridDefaultMpi());
	LatticePGG pgg(grid);
	LatticePGG avg(grid);
	avg = 0.;
  //
	std::string dir = "/home/ydzhao/cuth/three_point_config/32D_ama";
	//
	std::vector<int> trajs; 
	std::map<int, std::string> files;
	get_trajs(dir, trajs, files);

	for(int traj: trajs) {

		cout << "reading from: " << files[traj] << endl;
		readScidac(pgg, files[traj]);
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

	// writeScidac(avg, dir + "/average");
	
	// readScidac(avg, dir + "/average");
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
