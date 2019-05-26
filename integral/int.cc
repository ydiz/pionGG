#include <math.h>
#include <gsl/gsl_integration.h>
#include <cuba.h>
#include <string>
#include <vector>
#include <iostream>
#include <iomanip>

#include "../constants_macro.h"
#include "cuba_wrapper.h"



double f1(double p, void *params) {
  std::vector<double> paras = *(std::vector<double> *)params;
  double w = paras[0];
  double w0 = paras[1];
  double sin_pw, cos_pw;
  sincos(p * w, &sin_pw, &cos_pw);

  // return std::exp(-p * w0) / (M_PION - 2 * p) * (cos_pw - sin_pw / (p * w));
  return - 0.5 * std::exp(-p * w0) * (cos_pw - sin_pw / (p * w)); // f should not contain the pole term if using qaws
}

// p.s. integrate_qawc will add a term (p - pole) to the denominator; DO NOT write this factor in the function
double f2(double p, void *params) {
  std::vector<double> paras = *(std::vector<double> *)params;
  double w = paras[0];
  double w0 = paras[1];

  double sin_pw, cos_pw;
  sincos(p * w, &sin_pw, &cos_pw);

  return std::exp(-p * w0) / (M_PION + 2 * p) * (cos_pw - sin_pw / (p * w));
}


template<class T>
void integrate_qags(T f, std::vector<double> params, double lower, double upper, double epsabs, double epsrel, gsl_integration_workspace *w, double &result, double &error) {
  gsl_function F;
  F.function = f;
  F.params = &params;
  gsl_integration_qags (&F, lower, upper, epsabs, epsrel, 1000, w, &result, &error);
}

// calculate principal part integral
template<class T>
void integrate_qawc(T f, std::vector<double> params, double lower, double upper, double pole, double epsabs, double epsrel, gsl_integration_workspace *w, double &result, double &error) {
  gsl_function F;
  F.function = f;
  F.params = &params;
  gsl_integration_qawc(&F, lower, upper, pole, epsabs, epsrel, 1000, w, &result, &error);
}


struct Func {

  double p_interval;
	double pe = std::sqrt(0.5 * M_PION * 0.5 * M_PION - M_L *M_L);
  std::vector<double> paras; // [w, w0]

  std::vector<double> operator()(const std::vector<double>& v) const;
};

std::vector<double> Func::operator()(const std::vector<double>& v) const{

  std::vector<double> ans(1);

	double p = v[0] * p_interval; // p: [0, p_interval]
	double c = 2. * v[1] - 1.; // cos: [-1, 1]
	double Epe = std::sqrt(p*p + (M_PION/2.)*(M_PION/2.) - 2 * p * pe * c);

  double w = paras[0];
  double w0 = paras[1];

  double sin_pw, cos_pw;
  sincos(p * w, &sin_pw, &cos_pw);

  ans[0] =  std::exp(-Epe * w0) / (Epe * (- M_PION*M_PION + 4. * pe * pe * c * c)) * (cos_pw - sin_pw / (p * w));

  ans[0] *= 2. * p_interval;

  return ans;
}

template<class T>
double integrate_CUBA(T func, double eps_rel)
{
  std::vector<double> integral, error, prob;
  int nregions, neval, fail;

	integrateCuhre(integral, error, prob, nregions, neval, fail, 3, 1, func, 0., eps_rel);
  return integral[0];

}


gsl_integration_workspace * workspace = gsl_integration_workspace_alloc (1000);

int main (void)
{
  
  double ret1, ret2, ret3, error;
  double lower = 0.00001, upper = 30, epsrel = 1e-4;

  double beta = std::sqrt(1 - 4*M_L*M_L / (M_PION*M_PION));
  double log_beta = std::log((1 + beta) / (1 - beta));
  // double beta_coeff = log_beta / beta;

  Func f3;
  f3.p_interval = upper;

  for(double w=0; w<=28; w+=1.) 
    for(int w0=0; w0<=16; ++w0) {
      double ret;
      std::vector<double> params {w, double(w0)};
      std::cout << "w=" << w << " w0=" << w0 << std::endl;

      f3.paras = params;

      if(w==0) ret =0.;
      else {
        integrate_qawc(f1, params, lower, upper, M_PION/2., 0, epsrel, workspace, ret1, error);
        integrate_qags(f2, params, lower, upper, 0, epsrel, workspace, ret2, error);
        ret3 = integrate_CUBA(f3, epsrel);
        // ret = - std::exp(0.5 * M_PION * w0) / w / w * M_PI / 4. / (f3.pe * M_PION) * beta_coeff * ret1 
        //       + std::exp(-0.5 * M_PION * w0) / w / w * M_PI / 4. / (f3.pe * M_PION) * beta_coeff * ret2
        //       + M_PI / 2.0 / w / w * ret3;
        ret = - std::exp(0.5 * M_PION * w0) / w / w * M_PI / (f3.pe * M_PION) * log_beta * ret1 
              + std::exp(-0.5 * M_PION * w0) / w / w * M_PI / (f3.pe * M_PION) * log_beta * ret2
              + 2.0 * M_PI / w / w * ret3;
      }
      std::cout << std::setprecision(10) << "integral = " << ret  << std::endl;

  }

  gsl_integration_workspace_free (workspace);

  return 0;
}
