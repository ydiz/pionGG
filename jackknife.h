#pragma once

#include <Grid/Grid.h>
#include "lep_para.h"

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

  std::string ensemble;
  double M_h; // mass of pion/kaon in lattice unit
  double N_h; // normalization of wall src operator. N_h = <0 | pi(0) | pi>
  double Z_V;
  double hadron_coeff;
  std::string target;
  int traj_start, traj_end, traj_sep, traj_num;
  int time_cutoff_start, time_cutoff_end, time_cutoff_num;
  Lep_para lep_para;

  void get_three_point(LatticePGG &three_point, int traj);
  std::vector<double> get_result_with_cutoff(const LatticePGG &hadronic, const LatticePGG &leptonic);
  
};


void Jack_para::get_three_point(LatticePGG &three_point, int traj) {
  if(ensemble == "Pion_32ID") {
    std::string file = three_point_exact_path(traj); 
    read_cheng_PGG(three_point, file); // read 
  }
  else if(ensemble == "Pion_32IDF") {
    std::string file = three_point_path_32IDF(traj);
    read_luchang_PGG(three_point, file); // FIXME: change this after cheng generated his three point functions
  }
  else assert(0);
}

std::vector<double> Jack_para::get_result_with_cutoff(const LatticePGG &three_point, const LatticePGG &leptonic) {
  if(target=="form_factor") return form_factor(three_point, leptonic);
  else if(target == "real" || target == "real_CUBA3d" || target=="imag_analytic" || target == "imag_CUBA3d") return calculate_decay_rate_cutoff(three_point, leptonic, lep_para.lep_coef(), hadron_coeff, M_h);
  else assert(0);
}








}}
