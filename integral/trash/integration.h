#pragma once

#include <cuba.h>
#include <iostream>
#include <vector>
#include <cstring>
#include <cassert>
#include <math.h>
#include "cuba_wrapper.h"
#include "constants_macro.h"

struct Int_para {
	double eps_rel;
	double p_interval;
	double term1_eps;
	std::string algorithm;
	std::string term;
	int x_for_p1;

  double M_hadron;
  double M_lepton;
};

const double PI = 3.1415926535897932384626433832795;


class Func {

	std::vector<double> (Func::*pf)(const std::vector<double> &) const;
	const double M_hadron; // mass of pion/kaon
	const double M_lepton; // mass of electron/muon
	double pe; // momentum of outgoing electron/muon
	double eps; // eps for calculating the principal part of term 1

public:
	double p_interval;
	std::vector<double> v; // v is a three-vector; v = (sqrt(x^2 + y^2), z, |t|)

	Func(const Int_para& para);

	std::vector<double> operator()(const std::vector<double>& vx) const { return (this->*pf)(vx);	}
	
	std::vector<double> integral_total(const std::vector<double>& vx) const;
	std::vector<double> integral_total_cos(const std::vector<double>& vx) const;
	std::vector<double> spherical_integral_term_1(const std::vector<double>& vx) const; 
	std::vector<double> spherical_integral_term_2(const std::vector<double>& vx) const; 
	std::vector<double> spherical_integral_term_3(const std::vector<double>& vx) const; 

	std::vector<double> integral_total_with_pe0(const std::vector<double>& vx) const;
	std::vector<double> integral_total_with_p3(const std::vector<double>& vx) const;
	std::vector<double> integral_total_with_p1(const std::vector<double>& vx) const;
};

Func::Func(const Int_para& para) : M_hadron(para.M_hadron), M_lepton(para.M_lepton) {

	p_interval = para.p_interval;
	eps = para.term1_eps;
	pe = sqrt(0.5 * M_hadron * 0.5 * M_hadron - M_lepton*M_lepton);

	// if(para.term == "term1") pf = &Func::spherical_integral_term_1;
	// else if(para.term == "term2") pf = &Func::spherical_integral_term_2;
	// else if(para.term == "term3") pf = &Func::spherical_integral_term_3;
  // if(para.term == "all") pf = &Func::integral_total;
	// else if(para.term == "all_with_pe0") pf = &Func::integral_total_with_pe0;
	if(para.term == "with_p3") pf = &Func::integral_total_with_p3;
	else if(para.term == "with_p1") pf = &Func::integral_total_with_p1;
	else assert(0);

}



std::vector<double> Func::integral_total_with_p3(const std::vector<double>& vx) const
{
  std::vector<double> ans(1);
	double p = vx[0] * p_interval; // p: [0, p_interval]
	double c = 2. * vx[1] - 1.; // cos: [-1, 1]
	double phi = vx[2] * 2. * PI; // phi: [0, 2PI]

	double term1, term2, term3;

	double tmp = sqrt(1 - c*c)*cos(phi)*v[0] + c*v[1]; // \hat{p} \dot \vec{w} // performace can improve ~25%
	// term 2
  term2 = 4 * PI * p_interval * p * c * exp(- p * v[2]) * sin( p * tmp ) / ((p + M_hadron * 0.5) * (M_hadron * 0.5 + pe * c));

	// term 3
	double Epe = sqrt(p*p + (M_hadron/2.)*(M_hadron/2.) - 2 * p * pe * c);
  term3 = 4 * PI * p_interval * p * c * exp(- Epe * v[2]) * sin( p * tmp ) / (Epe * (M_hadron * 0.5 + pe * c) * (-M_hadron * 0.5 + pe * c));

	// term 1
	p = vx[0] * (M_hadron*0.5 - eps); // p: [0, M_hadron/2 - eps]

  term1 = 4 * PI * (M_hadron*0.5 - eps) * p * c * exp(- p * v[2]) * sin( p * tmp ) / ((-p + M_hadron * 0.5) * (-M_hadron * 0.5 + pe * c));

	p = vx[0]*(p_interval - (M_hadron*0.5 + eps)) + (M_hadron*0.5 + eps);  // p: [M_hadron/2 + eps, p_interval]

	term1 += 4 * PI * (p_interval - (M_hadron*0.5 + eps)) * p * c * exp(- p * v[2]) * sin( p * tmp ) / ((-p + M_hadron * 0.5) * (-M_hadron * 0.5 + pe * c));

	// final result
	// ans[0] = term1 + term2 + term3;
	ans[0] = (1. / M_hadron) * exp(0.5*M_hadron*v[2]) * term1 + (1. / M_hadron) * exp(-0.5*M_hadron*v[2]) * term2 + term3;
// 
  return ans;
}



std::vector<double> Func::integral_total_with_p1(const std::vector<double>& vx) const
{
  std::vector<double> ans(1);
	double p = vx[0] * p_interval; // p: [0, p_interval]
	double c = 2. * vx[1] - 1.; // cos: [-1, 1]
	double sin_theta = sqrt(1 - c*c); // sin_theta
	double phi = vx[2] * 2. * PI; // phi: [0, 2PI]

	double term1, term2, term3;
	double sin_phi, cos_phi;
	sincos(phi, &sin_phi, &cos_phi);

	double tmp = sin_theta*cos_phi*v[0] + sin_theta*sin_phi*v[1] + c*v[2]; // \hat{p} \dot \vec{w} // performace can improve ~25%
	// term 2
  term2 = 4 * PI * p_interval * p * sin_theta * cos_phi * exp(- p * v[3]) * sin( p * tmp ) / ((p + M_hadron * 0.5) * (M_hadron * 0.5 + pe * c));

	// term 3
	double Epe = sqrt(p*p + (M_hadron/2.)*(M_hadron/2.) - 2 * p * pe * c);
  term3 = 4 * PI * p_interval * p * sin_theta * cos_phi * exp(- Epe * v[3]) * sin( p * tmp ) / (Epe * (M_hadron * 0.5 + pe * c) * (-M_hadron * 0.5 + pe * c));

	// term 1
	p = vx[0] * (M_hadron*0.5 - eps); // p: [0, M_hadron/2 - eps]

  term1 = 4 * PI * (M_hadron*0.5 - eps) * p * sin_theta * cos_phi * exp(- p * v[3]) * sin( p * tmp ) / ((-p + M_hadron * 0.5) * (-M_hadron * 0.5 + pe * c));

	p = vx[0]*(p_interval - (M_hadron*0.5 + eps)) + (M_hadron*0.5 + eps);  // p: [M_hadron/2 + eps, p_interval]

	term1 += 4 * PI * (p_interval - (M_hadron*0.5 + eps)) * p * sin_theta * cos_phi * exp(- p * v[3]) * sin( p * tmp ) / ((-p + M_hadron * 0.5) * (-M_hadron * 0.5 + pe * c));

	// final result
	// ans[0] = term1 + term2 + term3;
	ans[0] = (1. / M_hadron) * exp(0.5*M_hadron*v[3]) * term1 + (1. / M_hadron) * exp(-0.5*M_hadron*v[3]) * term2 + term3;
// 
  return ans;
}



void integrate(const Int_para &para, const Func &func)
{
  std::vector<double> integral, error, prob;
  int nregions, neval, fail;

	if(para.algorithm=="Cuhre") integrateCuhre(integral, error, prob, nregions, neval, fail, 3, 1, func, 0., para.eps_rel);
	else if(para.algorithm=="Divonne") integrateDivonne(integral, error, prob, nregions, neval, fail, 3, 1, func, 0., para.eps_rel);
	else assert(0); // wrong algorithm name

}


void run_p1(const Int_para &para) {

	Func func(para);

	int x = para.x_for_p1;
	int y_max = SPACE_LIMIT;
	int z_max = SPACE_LIMIT;
	int t_max = SPACE_LIMIT; 
	std::cout << "x: " << x << std::endl; 
	std::cout << "y: [" << 0 << ", " << y_max << "]"<< std::endl; 
	std::cout << "z: [" << 0  << ", " << z_max << "]"<< std::endl; 
	std::cout << "t: [" << 0 << ", " << t_max << "]"<< std::endl; 
	std::cout << std::string(20, '*') << std::endl;

  // p1 integral is symmetric with relative to y and z;  and anti-symmetric in x
	// for(int x=x_min; x<=x_max; ++x) 
	for(int y=0; y<=y_max; ++y) 
		for(int z=0; z<=z_max; ++z)
			for(int t=0; t<=t_max; ++t) { // t [0, L/4] // t appears only as absolute value
				func.v = {double(x), double(y), double(z), double(t)};
				std::cout << "v: " << "[" << func.v[0] << " " << func.v[1] << " " << func.v[2] << " " << func.v[3] << "]" << std::endl;

        if(x==0) { // when x=0, integral is 0 and cuba cannote calculate it
          std::cout << cur_time << "s	integral = " << 0. << std::endl;
          continue;
        }	

				if(t==0 || t==1){ // Cuba cannot calculate when t=0; t=1 is very slow to calculate
					std::cout << cur_time << "s	Cannot calculate t=0,1; I am actually calculating with t=4" << std::endl;
					func.v[3] = 4;
				}	
        if(t==2 || t==3){ // Cuba cannot calculate when t=0; t=1 is very slow to calculate
					std::cout << cur_time << "calculating t=2,3 is slow; I am actually calculating with t=4" << std::endl;
					func.v[3] = 4;
				}
				integrate(para, func);
	}

	// func.v = {1., 0., -1., 10.};
	// std::cout << "v: " << "[" << func.v[0] << " " << func.v[1] << " " << func.v[2] << " " << func.v[3] << "]" << std::endl;
	// integrate(para, func);
}


void run_p3(const Int_para &para) {
	Func func(para);

	// func.v = {0, 0, 1};
	// integrate(para, func);
	// func.v = {sqrt(2), 1, 1};
  // return;

	int r_max = int(SPACE_LIMIT*std::sqrt(2)) + 1, z_max = SPACE_LIMIT, t_max = TIME_LIMIT;
	std::cout << "r: [" << 0 << ", " << r_max << "]"<< std::endl; 
	std::cout << "z: [" << 0 << ", " << z_max << "]"<< std::endl; 
	std::cout << "t: [" << 0 << ", " << t_max << "]"<< std::endl; 
	std::cout << std::string(20, '*') << std::endl;

	for(int r=0; r<=r_max; ++r) // r is distance, always positive
		for(int z=0; z<=z_max; ++z)
			for(int t=0; t<=t_max; ++t) { // t [0, L/4] // t appears only as absolute value

				func.v = {double(r), double(z), double(t)};
				std::cout << "v: " << "[" << func.v[0] << " " << func.v[1] << " " << func.v[2] << "]" << std::endl;

        // for sin, when z=0, integral is 0
        if(z==0) {
          std::cout << cur_time << "s	integral = " << 0. << std::endl;
          continue;
        }	

				if(t == 0 || t==1 || t==2){
          std::cout << cur_time << "s	Cannot calculate t=0; when r and z are large, t=1 is very slow to calculate; t=2 is not hard to do, but t's value does not affect the integral very much. So I am actually calculating with t=3" << std::endl;
          func.v[2] = 3;
        }

        integrate(para, func);
  }
}
