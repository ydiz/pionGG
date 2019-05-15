#include <Grid/Grid.h>
#include <dirent.h>
#include "constants_macro.h"
#include "io.h"
#include "pGG.h"
#include "read_leptonic.h"
#include "amplitude.h"

using namespace std;
using namespace Grid;
using namespace Grid::QCD;

const std::vector<int> gcoor({32, 32, 32, 64});

int main(int argc, char* argv[])
{

  Grid_init(&argc, &argv);

  GridCartesian * grid = SpaceTimeGrid::makeFourDimGrid(gcoor, GridDefaultSimd(Nd,vComplex::Nsimd()), GridDefaultMpi());

  int space_cutoff = 16;
  int time_cutoff = 16;

	LatticePGG leptonic(grid);

	LatticePGG avg(grid);
	// readScidac(avg, "./lat_config/average_three_point_sloppy");
	// readScidac(avg, "./lat_config/hadronic_wall_on_right/average_three_point_sloppy");
	readScidac(avg, "./lat_config/average_three_point_exact");

	// std::string filename_p3 = "/home/ydzhao/cuth/cuba_integration/results/p30e10-5Cuhre_with_p3/data.txt";
	// std::string filename_p1 = "/home/ydzhao/cuth/cuba_integration/results/p30e10-4Cuhre_with_p1/data.txt";
  std::string filename_p3 = "/projects/HadronicLight_4/yidizhao/cooley/pionGG/integrals/p30e10-5Cuhre_with_p3/data.txt";
  std::string filename_p1 = "/projects/HadronicLight_4/yidizhao/cooley/pionGG/integrals/p30e10-4Cuhre_with_p1/data.txt";
	get_leptonic(filename_p1, filename_p3, leptonic, SPACE_LIMIT, TIME_LIMIT);

  calculate_decay_rate(avg, leptonic, space_cutoff, time_cutoff); 


	Grid_finalize();
  return 0;
}
