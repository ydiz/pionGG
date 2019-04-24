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

int get_t(const std::string &path) {
	return std::stoi(path.substr(path.find("=") + 1));
}

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


	// int traj_start = 680, traj_num = 70;
	// std::vector<int> trajs(traj_num);
	std::vector<int> trajs {1370};
	// for(int i=0; i<trajs.size(); ++i) trajs[i] = traj_start + i * 10;

	cout << "trajs: " << endl;
	cout << trajs << endl;

	// std::vector<double> average_corr(gcoor[3], 0.);

	for(int traj: trajs) {

		std::string gauge_transform_path = gauge_transform_path_32D(traj);
		cout << "Load Gauge Transform And Get Inv: " <<  gauge_transform_path << endl;
		assert(dirExists(gauge_transform_path));
		GaugeTransform gtinv;
		{
			GaugeTransform gt;
			dist_read_field(gt, gauge_transform_path);
			to_from_big_endian_64(get_data(gt)); //FIXME: why? try commenting it
			gt_inverse(gtinv, gt);
		}
		std::vector<int> coor{1,13,23,47};
		print_qlat_field_site(gtinv, coor);


		std::vector<int> ts;
		std::map<int, std::string> wall_subdirs;
		std::string wall_src_path = wall_path_32D_sloppy(traj);
		get_ts(wall_src_path, ts, wall_subdirs);
		
		std::vector<std::vector<int>> xgs;
		std::map<std::vector<int>, std::string> point_subdirs;
		std::string point_src_path = point_path_32D(traj);
		get_xgs(point_src_path, xgs, point_subdirs);

		xgs.clear(); xgs.push_back({0, 4, 2, 6});

		for(const auto &xg: xgs) {

			// assert(dirExists(subdirs[xg]));
			// assert(read_mpi_coor(subdirs[xg]) == mpi_coor);
			cout << "xg of point src: " << xg << endl;
			cout << "directory name: " << point_subdirs[xg] << endl;

			qlat::Propagator4d qlat_point_prop;
			dist_read_field_double_from_float(qlat_point_prop, point_subdirs[xg]);
			std::cout << "Finished reading point propagator!" << std::endl;

			int wall_t = xg[3] + 10; 
			cout << "t of wall src: " << wall_t << endl;
			cout << "directory name: " << wall_subdirs[wall_t] << endl;

			qlat::Propagator4d qlat_wall_prop;
			dist_read_field_double_from_float(qlat_wall_prop, wall_subdirs[wall_t]);
			std::cout << "Finished reading wall propagator!" << std::endl;

			std::vector<int> coor{1,13,23,47};
			print_qlat_field_site(qlat_point_prop, coor);
			print_qlat_field_site(qlat_wall_prop, coor);

			// break;

			// grid_convert(prop, qlat_prop);
			// // print_grid_field_site(prop, coor);
      //
			// std::vector<typename LatticePropagator::vector_object::scalar_object> slice_sum;
			// sliceSum(prop, slice_sum, Tdir);
      //
			// for(int i=0; i<slice_sum.size(); ++i) {
			// 	int sep = (i < t) ? i - t + gcoor[3] : i - t;
			// 	double tmp = TensorRemove(trace(slice_sum[i] * adj(slice_sum[i]))).real();
			// 	corr[sep] += tmp;
			// }

		}

		// for(auto &x: corr) x /= gcoor[3];
		// cout << "Wall to Wall Correlator[traj=" << traj << "]:"<< endl;
		// cout << corr << endl;
    //
		// for(int i=0; i<corr.size(); ++i) average_corr[i] += corr[i];
	}

	// for(auto &x: average_corr) x /= trajs.size();
	// cout << std::string('*', 30) << endl;
	// cout << "average wall to wall correlator over "<< trajs.size() << " trajectoies:" << endl;
	// cout << average_corr << endl;

  end();

  return 0;
}
