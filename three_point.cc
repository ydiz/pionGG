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

  GridCartesian *grid = SpaceTimeGrid::makeFourDimGrid(gcoor, GridDefaultSimd(Nd,vComplex::Nsimd()), mpi_coor); 

	// int traj_start = 680, traj_num = 70;
	// std::vector<int> trajs(traj_num);
	// for(int i=0; i<trajs.size(); ++i) trajs[i] = traj_start + i * 10;
	std::vector<int> trajs {1370};

	cout << "trajs: " << endl;
	cout << trajs << endl;

	Gamma gamma5(Gamma::Algebra::Gamma5);
	// Gamma::gmu[0], Gamma::gmu[1], Gamma::gmu[0], Gamma::gmu[1];// GammaX, GammaY, GammaZ, GammaT

	for(int traj: trajs) {

		std::string gauge_transform_path = gauge_transform_path_32D(traj);

		std::vector<int> ts;
		std::map<int, std::string> wall_subdirs;
		std::string wall_src_path = wall_path_32D_sloppy(traj);
		get_ts(wall_src_path, ts, wall_subdirs);
		
		std::vector<std::vector<int>> xgs;
		std::map<std::vector<int>, std::string> point_subdirs;
		std::string point_src_path = point_path_32D(traj);
		get_xgs(point_src_path, xgs, point_subdirs);

    // read gauge transformation
		cout << "Load Gauge Transform And Get Inv: " <<  gauge_transform_path << endl;
		assert(dirExists(gauge_transform_path));
		GaugeTransform qlat_gtinv;
		{
			GaugeTransform qlat_gt;
			dist_read_field(qlat_gt, gauge_transform_path);
			to_from_big_endian_64(get_data(qlat_gt)); 
			gt_inverse(qlat_gtinv, qlat_gt);
		}
		LatticeColourMatrix gt(grid);
		grid_convert(gt, qlat_gtinv);
		// std::vector<int> coor{1,13,23,47};
		// print_grid_field_site(gt, coor);

    // read wall src
		std::vector<LatticePropagator> wall_props(gcoor[Tdir], grid);
		for(int t=0; t<gcoor[Tdir]; ++t) {
			read_propagator(wall_props[t], wall_subdirs[t]);
			wall_props[t] = gt * wall_props[t];
		}
		// std::cout << "Finished reading wall propagator!" << std::endl;

		xgs.clear(); xgs.push_back({0, 4, 2, 6});

		for(const auto &xg: xgs) {

			// assert(dirExists(subdirs[xg]));
			// assert(read_mpi_coor(subdirs[xg]) == mpi_coor);
			cout << "xg of point src: " << xg << endl;
			cout << "directory name: " << point_subdirs[xg] << endl;

			LatticePropagator point_prop(grid);
			read_propagator(point_prop, point_subdirs[xg]);
			std::cout << "Finished reading point propagator!" << std::endl;

      // !! global peekSite, not local peekLocalSite
      // !! cannot use global peekSite inside parallel_for; every node must peek the same site
      std::vector<typename LatticePropagator::vector_object::scalar_object> wall_props_to_xg(gcoor[Tdir]);
      for(int i=0; i<gcoor[Tdir]; ++i)	peekSite(wall_props_to_xg[i], wall_props[i], xg); 

			LatticePGG ret(grid);
			parallel_for(int ss=0; ss<ret._grid->lSites(); ss++){
				std::vector<int> lcoor(4);
				ret._grid->LocalIndexToLocalCoor(ss, lcoor);

				std::vector<int> gcoor(4);
				std::vector<int> processor_coor;
				ret._grid->ProcessorCoorFromRank(ret._grid->ThisRank(), processor_coor);
				ret._grid->ProcessorCoorLocalCoorToGlobalCoor(processor_coor, lcoor, gcoor);
				//std::cout << processor_coor << lcoor << gcoor << std::endl;

				std::vector<int> xp = gcoor;
				std::vector<int> x = xg;

        // cheng's tsep

        int t_min = 10;
        int t_wall;
        int t_sep;
      
        // using namespace qlat;
        int diff = smod(x[3] - xp[3], 64);
        if (diff >= 0)
        {
          t_wall = qlat::mod(x[3] + t_min, 64);
        } else {
          t_wall = qlat::mod(xp[3] + t_min, 64);
        }
        t_sep = qlat::mod(t_wall - x[3], 64);

				typename LatticePropagator::vector_object::scalar_object wall_to_x, x_to_wall, wall_to_xp, xp_to_wall, x_to_xp, xp_to_x;
				typename LatticePGG::vector_object::scalar_object ret_site;

				// peekSite(wall_to_x, wall_props[t_wall], x); //global, not local
        wall_to_x = wall_props_to_xg[t_wall];
				x_to_wall = gamma5 * adj(wall_to_x) * gamma5;
				peekLocalSite(wall_to_xp, wall_props[t_wall], lcoor);
				xp_to_wall = gamma5 * adj(wall_to_xp) * gamma5;
				peekLocalSite(x_to_xp, point_prop, lcoor);
				xp_to_x = gamma5 * adj(x_to_xp) * gamma5;

				ret_site = 0.;
				for(int mu=0; mu<4; ++mu)
					for(int nu=0; nu<4; ++nu) {
						ret_site()()(mu, nu) = trace(Gamma::gmu[mu] * xp_to_x * Gamma::gmu[nu] * wall_to_xp * gamma5 * x_to_wall); // cheng's order
						// ret_site()()(mu, nu) = trace(Gamma::gmu[mu] * xp_to_x * Gamma::gmu[nu] * wall_to_xp * gamma5 * x_to_wall) + trace(Gamma::gmu[nu] * x_to_xp * Gamma::gmu[mu] * wall_to_x * gamma5 * xp_to_wall); // cheng's order
				}

				ret_site = ret_site * (1 /  std::exp(- M_PION * t_sep));

				pokeLocalSite(ret_site, ret, lcoor);

			}

			for(int mu=0; mu<4; ++mu) ret = Cshift(ret, mu, xg[mu]);

      ret = real(ret);
      ret = 2.0 * ret; // two diagrams: clockwise and anti-clockwise; they are conjugate to each other.

      writeScidac(ret, "xg=(0,4,2,6)");
			cout << ret << endl;
		} // end of xg loop

	}


  end();

  return 0;
}
