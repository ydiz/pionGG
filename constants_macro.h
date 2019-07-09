#pragma once

/////////////////////////////////////


std::string loop_path_24ID(int traj) {
  std::string path ="/home/ljin/application/Public/Muon-GM2-cc/jobs/convert/data/24D/light-minus-heavy-1/traj=0" + std::to_string(traj); 
  assert(dirExists(path));
  return path;
}



//////////////////////////////////////////////////




std::string point_path_strange_24ID(int traj) {
  std::string path = "/home/ljin/application/Public/Muon-GM2-cc/jobs/32D/discon-1/results/prop-hvp ; results=" + std::to_string(traj) + "/huge-data/prop-point-src";
  assert(dirExists(path));
  return path;
}

std::string wall_path_ud_48I(int traj, int t) {
   std::string path = "/home/ljin/application/Public/Qlat-CPS-cc/jobs/48I/wall-src/results/48I-0.00078/results=" + std::to_string(traj) + "/huge-data/wall_src_propagator/t=" + std::to_string(t);
  assert(dirExists(path));
  return path;
}

std::string wall_path_strange_32D(int traj, int t) {
  std::string path = "/home/ljin/application/Public/Qlat-CPS-cc/jobs/wall-src-strange/results/32D-0.00107/results="+ std::to_string(traj) + "/huge-data/wall_src_propagator/strange ; t=" + std::to_string(t);
  assert(dirExists(path));
  return path;
}

// the same as wall_path_32D_sloppy
std::string wall_path_ud_32ID(int traj, int t) {
  std::string path = "/home/ljin/application/Public/Qlat-CPS-cc/jobs/32D/wall-src/results/32D-0.00107/results=" + std::to_string(traj) + "/huge-data/wall_src_propagator/t=" + std::to_string(t);
  assert(dirExists(path));
  return path;
}

std::string wall_path_24ID(int traj, int t) {
  std::string path = "/home/ljin/application/Public/Qlat-CPS-cc/jobs/24D/wall-src/results/results=" + std::to_string(traj) + "/huge-data/wall_src_propagator/t=" + std::to_string(t);
  assert(dirExists(path));
  return path;
}

std::string three_point_24ID(int traj) {
  std::string path = "/gpfs/mira-fs1/projects/HadronicLight_4/ctu/hlbl/hlbl-pion/TwoPointWallCorrField/24D-0.00107/ama/t-min=0010/results=" + std::to_string(traj) + "/avg ; type=0";    
  assert(dirExists(path));
  return path;
}


std::string three_point_32ID(int traj) {
	std::string str_traj = std::to_string(traj);
	if(str_traj.size()==3) str_traj = "0" + str_traj;

  std::string path = "/gpfs/mira-fs1/projects/HadronicLight_4/ctu/hlbl/hlbl-pion/TwoPointWallCorrField/32D-0.00107/ama/t-min=0010/results=" + str_traj + "/avg ; type=0";    
  assert(dirExists(path));
  return path;
}


std::string gauge_transform_path_32D(int traj) {
	return "/home/ljin/application/Public/Qlat-CPS-cc/jobs/32D/wall-src/results/32D-0.00107/results=" + std::to_string(traj) + "/huge-data/gauge-transform";
}

std::string gauge_transform_path_24D(int traj) {
	return "/home/ljin/application/Public/Qlat-CPS-cc/jobs/24D/wall-src/results/results=" + std::to_string(traj) + "/huge-data/gauge-transform";
}


// std::string wall_path_32DF(int traj) {
//   return "/home/ljin/application/Public/Qlat-CPS-cc/jobs/32Dfine/wall-src/results/32Dfine-0.0001/results=" + std::to_string(traj) + "/huge-data/wall_src_propagator";
// }
//
// std::string wall_path_32D_sloppy(int traj) {
// 	return "/home/ljin/application/Public/Qlat-CPS-cc/jobs/32D/wall-src/results/32D-0.00107/results=" + std::to_string(traj) + "/huge-data/wall_src_propagator";
// }
//
// std::string wall_path_32D_exact(int traj) {
// 	return "/home/ljin/application/Public/Qlat-CPS-cc/jobs/32D/wall-src-exact-2/results/32D-0.00107/results="+ std::to_string(traj) + "/huge-data/wall_src_propagator";
// }





std::string point_path_32D(int traj) {
	return "/home/ljin/application/Public/Muon-GM2-cc/jobs/32D/discon-1/results/prop-hvp ; results=" + std::to_string(traj) + "/huge-data/prop-point-src";
}

std::string gauge_transform_path(int traj) {
	return "/home/ljin/application/Public/Qlat-CPS-cc/jobs/32D/wall-src/results/32D-0.00107/results=" + std::to_string(traj) + "/huge-data/gauge-transform";
}

// std::string three_point_sloppy_path(int traj) {
// 	std::string str_traj = std::to_string(traj);
// 	if(str_traj.size()==3) str_traj = "0" + str_traj;
//
// 	std::string path = "/gpfs/mira-fs1/projects/HadronicLight_4/ctu/hlbl/hlbl-pion/ThreePointCorrField/32D-0.00107/sloppy/results=" + str_traj + "/t-min=0010/avg ; type=0 ; accuracy=0";
// 	
// 	return path;
// }
//
// std::string three_point_exact_path(int traj) {
// 	std::string str_traj = std::to_string(traj);
// 	if(str_traj.size()==3) str_traj = "0" + str_traj;
//
// 	std::string path = "/gpfs/mira-fs1/projects/HadronicLight_4/ctu/hlbl/hlbl-pion/ThreePointCorrField/32D-0.00107/ama/results=" + str_traj + "/ama ; t-min=0010/avg ; type=0";
// 	
// 	return path;
// }



std::string three_point_path_32IDF(int traj) {

  // std::string path = "/home/ljin/application/Public/Qlat-CPS-cc/jobs/em-corr/results/32Dfine-0.0001/results="+ std::to_string(traj) + "/contraction-with-point/pion_gg/decay_cheng";
  std::string path = "/home/ljin/application/Public/Qlat-CPS-cc/jobs/em-corr/results/32Dfine-0.0001/results="+ std::to_string(traj) + "/contraction-with-point-old/pion_gg/decay_cheng";

  // std::cout << path << std::endl;
  assert(dirExists(path));
	
	return path;
}

std::string three_point_48I(int traj) {

  std::string path = "/home/ljin/application/Public/Qlat-CPS-cc/jobs/em-corr/results/48I-0.00078/results=" + std::to_string(traj) + "/contraction-with-point/pion_gg/decay_cheng";
  std::cout << path << std::endl;
  assert(dirExists(path));
	
	return path;
}


