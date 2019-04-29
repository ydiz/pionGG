#include <Grid/Grid.h>
#include <dirent.h>
#include "constants_macro.h"
#include "io.h"
#include "pGG.h"
#include "qlat_wrapper/qlat_wrapper.h"

using namespace std;
using namespace qlat;
using namespace Grid;
using namespace Grid::QCD;

const std::vector<int> gcoor({32, 32, 32, 64});

const std::vector<int> mpi_coor {1,1,1,8};
int main(int argc, char* argv[])
{

  begin(&argc, &argv, Coordinate(mpi_coor[0], mpi_coor[1], mpi_coor[2], mpi_coor[3])); // begin is defined in qlat/mpi.h
  Grid_init(&argc, &argv);

  GridCartesian * grid = SpaceTimeGrid::makeFourDimGrid(gcoor, GridDefaultSimd(Nd,vComplex::Nsimd()), mpi_coor);
	// cout << hadronic << endl;
  LatticePGG pgg(grid);
	string path = "/gpfs/mira-fs1/projects/HadronicLight_4/ctu/hlbl/hlbl-pion/ThreePointCorrField/32D-0.00107/sloppy/results=1370/t-min=0010/xg=(0,4,2,6) ; type=0 ; accuracy=0";
		qlat::PionGGElemField qlat_pgg;
		dist_read_field(qlat_pgg, path);
		cout << "Finished reading" << endl;
		grid_convert(pgg, qlat_pgg);

		cout << pgg << endl;


	Grid_finalize();
  return 0;
}
