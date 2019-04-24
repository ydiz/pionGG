#include <qlat/qlat.h>
#include <qlat/grid.h>
// #include "cheng.h"
// #include "util.h"

#include "qlat_wrapper/qlat_wrapper.h"

#include "constants_macro.h"

#include <dirent.h>
#include <sys/stat.h>

using namespace qlat;
using namespace std;
using namespace Grid;
using namespace Grid::QCD;

// std::string wall_src_path = "/home/ydzhao/cuth/luchang_propagators/32D_wall_src/results=1000/huge-data/wall_src_propagator";

int get_t(const std::string &path) {
	return std::stoi(path.substr(path.find("=") + 1));
}

// void get_trajs(const std::string &path, std::vector<int> &trajs, std::map<int, std::string> &subdirs) {
//     DIR *dir;
//     dir = opendir(path.c_str());
//     assert(dir!=NULL); // make sure directory exists
//     struct dirent *entry;
//
//     std::string filename;
//     while ((entry = readdir (dir)) != NULL) {
//       // printf ("%s\n", entry->d_name);
//       filename = std::string(entry->d_name);
// 	  if(filename.substr(0, 7) == "results") {
// 		  int traj = std::stoi(filename.substr( filename.find("=") + 1 )); // "results=1000";
// 		  trajs.push_back(traj);
// 		  subdirs.insert(std::pair<int, std::string>(traj, path + "/" + filename));
// 	  }
//     }
//     closedir (dir);
//
//     std::sort(trajs.begin(), trajs.end());
// }

void get_ts(const std::string &path, std::vector<int> &ts, std::map<int, std::string> &subdirs) {

		ts.clear();
		subdirs.clear();
    DIR *dir;
    dir = opendir(path.c_str());
    assert(dir!=NULL); // make sure directory exists
    struct dirent *entry;

    std::string subdir_name;
    while ((entry = readdir (dir)) != NULL) {
      // printf ("%s\n", entry->d_name);
      subdir_name = std::string(entry->d_name);
			if(subdir_name.substr(0, 2) == "t=") {
				int t = get_t(subdir_name); 
				ts.push_back(t);
				subdirs.insert(std::pair<int, std::string>(t, path + "/" + subdir_name));
			}
    }
    closedir (dir);

    std::sort(ts.begin(), ts.end());
}


std::vector<int> gcoor({32, 32, 32, 64});
std::vector<int> mpi_coor({1, 1, 1, 8});

int main(int argc, char* argv[])
{
  begin(&argc, &argv, Coordinate(mpi_coor[0], mpi_coor[1], mpi_coor[2], mpi_coor[3])); // begin is defined in qlat/mpi.h
  Grid_init(&argc, &argv);

  GridCartesian * grid = SpaceTimeGrid::makeFourDimGrid(gcoor, GridDefaultSimd(Nd,vComplex::Nsimd()), mpi_coor); 
	LatticePropagator prop(grid);


	int traj_start = 680, traj_num = 70;
	std::vector<int> trajs(traj_num);
	for(int i=0; i<trajs.size(); ++i) trajs[i] = traj_start + i * 10;

	cout << "trajs: " << endl;
	cout << trajs << endl;

	std::vector<double> average_corr(gcoor[3], 0.);

	for(int traj: trajs) {

		std::vector<int> ts;
		std::map<int, std::string> subdirs;
		std::string wall_src_path = wall_path_32D_sloppy(traj);
		get_ts(wall_src_path, ts, subdirs);

		std::vector<double> corr(gcoor[3], 0.);

		for(const auto &t: ts) {
			// cout << "t of wall src: " << t << endl;
			// cout << "directory name: " << subdirs[t] << endl;

			assert(dirExists(subdirs[t]));
			assert(read_mpi_coor(subdirs[t]) == mpi_coor);

			qlat::Propagator4d qlat_prop;
			dist_read_field_double_from_float(qlat_prop, subdirs[t]);
			// std::cout << "Finished reading propagator!" << std::endl;
			// std::vector<int> coor{1,13,23,47};
			// print_qlat_field_site(qlat_prop, coor);

			grid_convert(prop, qlat_prop);
			// print_grid_field_site(prop, coor);

			std::vector<typename LatticePropagator::vector_object::scalar_object> slice_sum;
			sliceSum(prop, slice_sum, Tdir);

			for(int i=0; i<slice_sum.size(); ++i) {
				int sep = (i < t) ? i - t + gcoor[3] : i - t;
				double tmp = TensorRemove(trace(slice_sum[i] * adj(slice_sum[i]))).real();
				corr[sep] += tmp;
			}

		}

		for(auto &x: corr) x /= gcoor[3];
		cout << "Wall to Wall Correlator[traj=" << traj << "]:"<< endl;
		cout << corr << endl;

		for(int i=0; i<corr.size(); ++i) average_corr[i] += corr[i];
	}

	for(auto &x: average_corr) x /= trajs.size();
	cout << std::string('*', 30) << endl;
	cout << "average wall to wall correlator over "<< trajs.size() << " trajectoies:" << endl;
	cout << average_corr << endl;

  end();

  return 0;
}
