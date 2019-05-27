// g++ -I/home/ydzhao/cuth/install/boost/include -L/home/ydzhao/cuth/install/boost/lib GF_para.cc -lboost_program_options
#pragma once

#include <stdlib.h>
#include <boost/program_options.hpp>
#include "jackknife.h"

namespace po = boost::program_options;

namespace Grid {
namespace QCD {

void init_para(int argc, char **argv, Jack_para &para)
{
  po::options_description desc("jackknife options");
  desc.add_options()("help", "help message")
                    ("ensemble", po::value<std::string>(&para.ensemble))
                    ("traj_start", po::value<int>(&para.traj_start)->default_value(1000))
                    ("traj_end", po::value<int>(&para.traj_end)->default_value(1000))
                    ("traj_sep", po::value<int>(&para.traj_sep)->default_value(10))
                    ("time_cutoff_start", po::value<int>(&para.time_cutoff_start)->default_value(1))
                    ("time_cutoff_end", po::value<int>(&para.time_cutoff_end)->default_value(16))
                    ("target", po::value<std::string>(&para.lep_para.target)->default_value(""))
                    ("file_p3", po::value<std::string>(&para.lep_para.file_p3)->default_value(""))
                    ("file_p1", po::value<std::string>(&para.lep_para.file_p1)->default_value(""))
                    ;

  po::variables_map vm;
  po::store(po::command_line_parser(argc, argv).options(desc).allow_unregistered().run(), vm); // allow additional command line options // command line options have higher priority
  po::store(po::parse_config_file<char>("jack_init.ini", desc), vm);
  po::notify(vm);

  para.traj_num = (para.traj_end - para.traj_start) / para.traj_sep + 1;
  para.time_cutoff_num = para.time_cutoff_end - para.time_cutoff_start + 1;

  para.target = para.lep_para.target;

  if(para.ensemble == "Pion_32ID") {
    para.M_h = 0.139474;
    para.N_h = 52.089753;
    para.Z_V = 0.7260;
  }
  else if(para.ensemble == "Pion_32IDF") {
    para.M_h = 0.10025; // FIXME
    para.N_h = 67.232746; // FIXME
    para.Z_V = 0.6885; // FIXME
  }
  else assert(0);

  para.hadron_coeff = 1./ (3 * std::sqrt(2)) * para.Z_V * para.Z_V * 2. * para.M_h / para.N_h;


  if(vm.count("help")) {
    std::cout << desc << std::endl;
    exit(0);
  }

	std::cout << std::string(20, '*') << std::endl;
  std::cout << "M_h: " << para.M_h << std::endl;
  std::cout << "N_h: " << para.N_h << std::endl;
  std::cout << "Z_V: " << para.Z_V << std::endl;
	std::cout << std::string(20, '*') << std::endl;
	std::cout << "target: " << para.lep_para.target << std::endl;
	std::cout << "file_p3: " << para.lep_para.file_p3 << std::endl;
	std::cout << "file_p1: " << para.lep_para.file_p1 << std::endl;
	std::cout << std::string(20, '*') << std::endl;
	std::cout << "traj_start: " << para.traj_start << std::endl;
	std::cout << "traj_end: " << para.traj_end << std::endl;
	std::cout << "traj_sep: " << para.traj_sep << std::endl;
	std::cout << "traj_num: " << para.traj_num << std::endl;
	std::cout << std::string(20, '*') << std::endl;
	std::cout << "time_cutoff_start: " << para.time_cutoff_start << std::endl;
	std::cout << "time_cutoff_end: " << para.time_cutoff_end << std::endl;
	std::cout << "time_cutoff_num: " << para.time_cutoff_num << std::endl;
	std::cout << std::string(20, '*') << std::endl;

}



}}
