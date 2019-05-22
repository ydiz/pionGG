#include "lep.h"

using namespace std;
using namespace Grid;
using namespace Grid::QCD;

const std::vector<int> gcoor({32, 32, 32, 64});

int main(int argc, char* argv[])
{

  Grid_init(&argc, &argv);

  GridCartesian * grid = SpaceTimeGrid::makeFourDimGrid(gcoor, GridDefaultSimd(Nd,vComplex::Nsimd()), GridDefaultMpi());

	std::string filename = "./integral/data.txt"; // integral with p3 in the numerator
  
	// int space_limit = 16, time_limit = 16;// x,y,z [-space_limit, space_limit]; t [-time_limit, time_limit]
	int space_limit = LEPTONIC_SPACE_LIMIT, time_limit = LEPTONIC_TIME_LIMIT;

	LatticePGG ret(grid);
	get_leptonic_new(filename, ret, space_limit, time_limit);
	// writeScidac(ret, "./lat_config/lat_leptonic");
  cout << ret << endl;

	Grid_finalize();
  return 0;
}
