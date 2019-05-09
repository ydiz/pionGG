#include "read_leptonic.h"

using namespace std;
using namespace Grid;
using namespace Grid::QCD;

const std::vector<int> gcoor({32, 32, 32, 64});

int main(int argc, char* argv[])
{

  Grid_init(&argc, &argv);

  GridCartesian * grid = SpaceTimeGrid::makeFourDimGrid(gcoor, GridDefaultSimd(Nd,vComplex::Nsimd()), GridDefaultMpi());

	// std::string filename_p3 = "/home/yidizhao/cooley/pionGG/integrals/p30e10-5Cuhre_with_p3/data.txt"; // integral with p3 in the numerator
	std::string filename_p3 = "/home/ydzhao/cuth/cuba_integration/results/p30e10-5Cuhre_with_p3/data.txt"; // integral with p3 in the numerator
	std::string filename_p1 = "/home/ydzhao/cuth/cuba_integration/results/p30e10-4Cuhre_with_p1/data.txt"; // integral with p1 in the numerator
  
	// int space_limit = 8, time_limit = 16;// x,y,z [-space_limit, space_limit]; t [0, time_limit]
	int space_limit = SPACE_LIMIT, time_limit = TIME_LIMIT;// x,y,z [-space_limit, space_limit]; t [0, time_limit]

	LatticePGG ret(grid);
	get_leptonic(filename_p1, filename_p3, ret, space_limit, time_limit);
	// writeScidac(ret, "./lat_config/lat_leptonic");
  cout << ret << endl;

	Grid_finalize();
  return 0;
}
