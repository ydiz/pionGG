#pragma once

#include <Grid/Grid.h>
#include <dirent.h>
#include "constants_macro.h"
#include "io.h"
#include "pGG.h"
#include "utils.h"


namespace Grid {
namespace QCD {

double linear_interpolation(double x, double x_lower, double y_lower, double x_upper, double y_upper) {
	return ( (x - x_lower) * y_upper +  (x_upper - x) * y_lower ) / (x_upper - x_lower);
}

using Data = std::vector<std::vector<double>>;

// double get_p3_site(const std::vector<int> &gcoor, double data[13][17][17], const std::vector<int> &gdim) {
double get_integral_site(const std::vector<int> &gcoor, const Data &data, const std::vector<int> &gdim, int space_limit, int time_limit) {

  // FIXME: delete this
	// int x = qlat::smod(gcoor[0], gdim[0]);
	// int y = qlat::smod(gcoor[1], gdim[1]);
	// int z = qlat::smod(gcoor[2], gdim[2]);
	// int t = qlat::smod(gcoor[3], gdim[3]);

	if( (x <= space_limit && x >= -space_limit) && (y <= space_limit && y >= -space_limit) && (z <= space_limit && z >= -space_limit) && (t<=time_limit && t >= -time_limit)) {
		double r = std::sqrt(x*x + y*y + z*z); 	
    // if(r==0) return 0.; // do not know why without this the site will be nan
		int floor = int(r), ceiling = int(r) + 1;
		double ret = linear_interpolation(r, floor, data[floor][std::abs(t)], ceiling, data[ceiling][std::abs(t)]);

    // return ret / r; //FIXME: r is divided here
    return ret; //FIXME: r is divided here
	}
	else return 0.;
}


void read_integrals(const std::string &filename, int space_limit, int time_limit, Data &data) 
{

	std::vector<int> shape {int(space_limit * std::sqrt(3)) + 2, time_limit + 1};
	assert(shape[0]==29 && shape[1]==17);	
  
  data.resize(shape[0]);
	for(int i = 0; i < shape[0]; i++)
	{
			data[i].resize(shape[1]);
	}

	std::ifstream f(filename);
	
	std::cout << "before loop " << filename << std::endl;
	for(int i=0; i<shape[0]; ++i)
		for(int j=0; j<shape[1]; ++j)
				f >> data[i][j];

	assert(!f.eof()); // make sure there aren't too few numbers
	double tmp;
	f>>tmp;
	assert(f.eof()); // make sure there aren't too many numbers
	
	f.close();
	std::cout << "Done reading integral from " << filename << std::endl;
}


void get_leptonic_new(const std::string &filename, LatticePGG &lat, int space_limit, int time_limit) {
	
	Data data;
	read_integrals(filename, space_limit, time_limit, data);

	parallel_for(int ss=0; ss<lat._grid->lSites(); ss++){

    std::vector<int> lcoor, gcoor;
    localIndexToLocalGlobalCoor(lat._grid, ss, lcoor, gcoor);

    gcoor = my_smod(gcoor, lat._grid->_fdimensions);

		double val = get_integral_site(gcoor, data, lat._grid->_fdimensions, space_limit, time_limit); 

		typename LatticePGG::vector_object::scalar_object m;
		m = 0.;
		m()()(0, 1) = Complex(val * gcoor[Zdir], 0); 
		m()()(0, 2) = Complex(-val * gcoor[Ydir], 0); // Minus sign comes from spinor matrix
		m()()(1, 2) = Complex(val * gcoor[Xdir], 0); 
		m()()(1, 0) = - m()()(0, 1);
		m()()(2, 0) = - m()()(0, 2);
		m()()(2, 1) = - m()()(1, 2);

		pokeLocalSite(m, lat, lcoor);
	}
}

}}
