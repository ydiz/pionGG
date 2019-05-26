#pragma once

#include <Grid/Grid.h>
#include "lep_para.h"

// void calculate_jackknife(const std::vector<double> &jackknife_results){
//
//   int traj_num = jackknife_results.size();
//
//   double jackknife_avg = 0.;
//   for(double x: jackknife_results) jackknife_avg += x;
//   jackknife_avg /= double(traj_num);
//
//   double jackknife_error = 0.;
//   for(double x: jackknife_results) jackknife_error += (x - jackknife_avg) * (x - jackknife_avg);
//   jackknife_error = std::sqrt(jackknife_error * (double(traj_num) - 1.) / double(traj_num));
//
//   std::cout << "jackknife average: " << jackknife_avg << std::endl;
//   std::cout << "jackknife error: " << jackknife_error << std::endl;
//   std::cout << std::string(20, '*') << std::endl;
//
// }
//
namespace Grid{
namespace QCD{

RealD mean(const std::vector<RealD>& data)
{
    int N = data.size();
    RealD mean(0.0);
    for(int i=0; i<N; ++i){ mean += data[i]; }
    return mean/RealD(N);
}

RealD jack_mean(const std::vector<RealD>& data, int sample)
{
    int N = data.size();
    RealD mean(0.0);
    for(int i=0; i<N; ++i){ if(i != sample){ mean += data[i]; } }
    return mean/RealD(N-1);
}

RealD jack_std(const std::vector<RealD>& jacks, RealD mean)
{
    int N = jacks.size();
    RealD std(0.0);
    for(int i=0; i<N; ++i){ std += std::pow(jacks[i]-mean, 2.0); }
    return std::sqrt(RealD(N-1)/RealD(N)*std);
}

std::vector<RealD> jack_stats(const std::vector<RealD>& data)
{
  int N = data.size();
  std::vector<RealD> jack_samples(N);
  std::vector<RealD> jack_stats(2);

  jack_stats[0] = mean(data);
  for(int i=0; i<N; i++){ jack_samples[i] = jack_mean(data,i); }
  jack_stats[1] = jack_std(jack_samples, jack_stats[0]);

  std::cout << "jackknife average: " << jack_stats[0] << std::endl;
  std::cout << "jackknife error: " << jack_stats[1] << std::endl;
  std::cout << std::string(20, '*') << std::endl;

  return jack_stats;
}


struct Jack_para {

  std::string target;
  int traj_start, traj_end, traj_sep, traj_num;
  int time_cutoff_start, time_cutoff_end, time_cutoff_num;
  Lep_para lep_para;

  std::vector<double> get_result_with_cutoff(const LatticePGG &hadronic, const LatticePGG &leptonic);
  
};

std::vector<double> Jack_para::get_result_with_cutoff(const LatticePGG &three_point, const LatticePGG &leptonic) {
  if(target=="form_factor") return form_factor(three_point, leptonic);
  else if(target == "real" || target == "real_CUBA3d" || target=="imag_analytic" || target == "imag_CUBA3d") return calculate_decay_rate_cutoff(three_point, leptonic, lep_para.lep_coef());
  else assert(0);
}








}}
