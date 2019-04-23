#include <Grid/Grid.h>
#include <dirent.h>
#include "constants_macro.h"
#include "io.h"
#include "pGG.h"


#define PI_MASS 0.135

namespace Grid {
namespace QCD {

double linear_interpolation(double x, double x_lower, double y_lower, double x_upper, double y_upper) {
	return ( (x - x_lower) * y_upper +  (x_upper - x) * y_lower ) / (x_upper - x_lower);
}

using Data = std::vector<std::vector<std::vector<double>>>;

// double get_p3_site(const std::vector<int> &gcoor, double data[13][17][17], const std::vector<int> &gdim) {
double get_integral_site(const std::vector<int> &gcoor, const Data &data, const std::vector<int> &gdim, int space_limit, int time_limit) {

	int x = smod(gcoor[0], gdim[0]);
	int y = smod(gcoor[1], gdim[1]);
	int z = smod(gcoor[2], gdim[2]);
	int t = smod(gcoor[3], gdim[3]);

	if( (x <= space_limit && x >= -space_limit) && (y <= space_limit && y >= -space_limit) && (z <= space_limit && z >= -space_limit) && (t<=time_limit && t >= -time_limit)) {
		double r = std::sqrt(x*x + y*y); 	
		int floor = int(r), ceiling = int(r) + 1;
		return linear_interpolation(r, floor, data[floor][z+space_limit][std::abs(t)], ceiling, data[ceiling][z+space_limit][std::abs(t)]);
	}
	else return 0.;
}

// double data[13][17][17]; //r[0, 12], z[-8, 8], t [0, 16]
void read_integrals(const std::string &filename, int space_limit, int time_limit, Data &data) 
{
	// x,y,z [-space_limit, space_limit] -> r [0, int(sqrt(x*x + y*y)) + 1]; t [0, time_limit]
	std::vector<int> shape {int(space_limit * std::sqrt(2)) + 2, 2 * space_limit + 1, time_limit + 1};
	// for(auto x: shape) std::cout << x << std::endl;
	// assert(shape[0]==13 && shape[1]==17 && shape[2]==17);	data.resize(shape[0]);
	assert(shape[0]==24 && shape[1]==33 && shape[2]==17);	data.resize(shape[0]);

	for(int i = 0; i < shape[0]; i++)
	{
			data[i].resize(shape[1]);
			for(int j = 0; j < shape[1]; j++)
					data[i][j].resize(shape[2]);
	}

	std::ifstream f(filename);
	
	for(int i=0; i<shape[0]; ++i)
		for(int j=0; j<shape[1]; ++j)
			for(int k=0; k<shape[2]; ++k) {
				f >> data[i][j][k];
			}

	assert(!f.eof()); // make sure there aren't too few numbers
	double tmp;
	f>>tmp;
	assert(f.eof()); // make sure there aren't too many numbers
	
	f.close();
	std::cout << "Done reading integral from " << filename << std::endl;
	// for(int i=0; i<shape[0]; ++i)
	// 	for(int j=0; j<shape[1]; ++j)
	// 		for(int k=0; k<shape[2]; ++k) {
	// 			std::cout << i << " "<<j <<" " <<k << ": " << data[i][j][k] << std::endl;
	// 		}
}

void get_leptonic(const std::string &filename_0, const std::string &filename_p3, LatticePGG &lat, int space_limit, int time_limit) {
	
	Data data_0, data_p3;
	read_integrals(filename_0, space_limit, time_limit, data_0);
	read_integrals(filename_p3, space_limit, time_limit, data_p3);

	parallel_for(int ss=0; ss<lat._grid->lSites(); ss++){
		std::vector<int> lcoor(4);
		lat._grid->LocalIndexToLocalCoor(ss, lcoor);

		std::vector<int> gcoor(4);
		std::vector<int> processor_coor;
		lat._grid->ProcessorCoorFromRank(lat._grid->ThisRank(), processor_coor);
		lat._grid->ProcessorCoorLocalCoorToGlobalCoor(processor_coor, lcoor, gcoor);
		//std::cout << processor_coor << lcoor << gcoor << std::endl;

		double val_0 = get_integral_site(gcoor, data_0, lat._grid->_fdimensions, space_limit, time_limit); 
		double val_p3 = get_integral_site(gcoor, data_p3, lat._grid->_fdimensions, space_limit, time_limit); 
		
		// typename LatticeComplex::vector_object::scalar_object m;
		typename LatticePGG::vector_object::scalar_object m;
		m = 0.;
		m()()(0, 1) = Complex(val_p3, 0); 
		m()()(2, 3) = Complex( - PI_MASS * 0.5 * val_0, 0);
		m()()(1, 0) = - m()()(0, 1);
		m()()(3, 2) = - m()()(2, 3);
		
		if(gcoor[Tdir] > lat._grid->_fdimensions[Tdir] / 2) m = -m; // FIXME: make sure this is right; when w-> -w, leptonic indeed becomes negative

		pokeLocalSite(m, lat, lcoor);
	}

	writeScidac(lat, "lat_leptonic");
	std::cout << lat << std::endl;
}

}}

using namespace std;
using namespace Grid;
using namespace Grid::QCD;

const std::vector<int> gcoor({32, 32, 32, 64});

int main(int argc, char* argv[])
{

  Grid_init(&argc, &argv);

  GridCartesian * grid = SpaceTimeGrid::makeFourDimGrid(gcoor, GridDefaultSimd(Nd,vComplex::Nsimd()), GridDefaultMpi());
	// LatticePGG avg(grid);
	// avg = 0.;

	std::string filename_0 = "/home/ydzhao/cuth/cuba_integration/results/p30e10-5Cuhre/data.txt"; // integral without p in the numerator
	std::string filename_p3 = "/home/ydzhao/cuth/cuba_integration/results/p30e10-5Cuhre_with_p3/data.txt"; // integral with p3 in the numerator

	LatticePGG ret(grid);
	// LatticeComplex ret(grid);

	// int space_limit = 8, time_limit = 16;// x,y,z [-space_limit, space_limit]; t [0, time_limit]
	int space_limit = SPACE_LIMIT, time_limit = TIME_LIMIT;// x,y,z [-space_limit, space_limit]; t [0, time_limit]
	get_leptonic(filename_0, filename_p3, ret, space_limit, time_limit);

	Grid_finalize();
  return 0;
}
