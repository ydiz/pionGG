#pragma once

#define LEPTONIC_TIME_LIMIT 16
#define LEPTONIC_SPACE_LIMIT 16

#define M_PION 0.139474 
#define M_L  0.000528 // calculate by 0.000511 * 0.139474 / 0.13497

std::string gauge_transform_path_32D(int traj) {
	return "/home/ljin/application/Public/Qlat-CPS-cc/jobs/32D/wall-src/results/32D-0.00107/results=" + std::to_string(traj) + "/huge-data/gauge-transform";
}

std::string gauge_transform_path_24D(int traj) {
	return "/home/ljin/application/Public/Qlat-CPS-cc/jobs/24D/wall-src/results/results=" + std::to_string(traj) + "/huge-data/gauge-transform";
}

std::string wall_path_32D_sloppy(int traj) {
	return "/home/ljin/application/Public/Qlat-CPS-cc/jobs/32D/wall-src/results/32D-0.00107/results=" + std::to_string(traj) + "/huge-data/wall_src_propagator";
}

std::string wall_path_32D_exact(int traj) {
	return "/home/ljin/application/Public/Qlat-CPS-cc/jobs/32D/wall-src-exact-2/results/32D-0.00107/results="+ std::to_string(traj) + "/huge-data/wall_src_propagator";
}

std::string wall_path_24D_sloppy(int traj) {
	return "/home/ljin/application/Public/Qlat-CPS-cc/jobs/24D/wall-src/results/results=" + std::to_string(traj) + "/huge-data/wall_src_propagator";
}

std::string point_path_32D(int traj) {
	return "/home/ljin/application/Public/Muon-GM2-cc/jobs/32D/discon-1/results/prop-hvp ; results=" + std::to_string(traj) + "/huge-data/prop-point-src";
}

// FIXME: ? do not distinguish between sloppy and exact
std::string gauge_transform_path(int traj) {
	return "/home/ljin/application/Public/Qlat-CPS-cc/jobs/32D/wall-src/results/32D-0.00107/results=" + std::to_string(traj) + "/huge-data/gauge-transform";
}

std::string three_point_sloppy_path(int traj) {
	std::string str_traj = std::to_string(traj);
	if(str_traj.size()==3) str_traj = "0" + str_traj;

	std::string path = "/gpfs/mira-fs1/projects/HadronicLight_4/ctu/hlbl/hlbl-pion/ThreePointCorrField/32D-0.00107/sloppy/results=" + str_traj + "/t-min=0010/avg ; type=0 ; accuracy=0";
	
	return path;
}

std::string three_point_exact_path(int traj) {
	std::string str_traj = std::to_string(traj);
	if(str_traj.size()==3) str_traj = "0" + str_traj;

	std::string path = "/gpfs/mira-fs1/projects/HadronicLight_4/ctu/hlbl/hlbl-pion/ThreePointCorrField/32D-0.00107/ama/results=" + str_traj + "/ama ; t-min=0010/avg ; type=0";
	
	return path;
}

