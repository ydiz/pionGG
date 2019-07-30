#pragma once

/////////////////////////////////////


// std::string loop_path_24ID(int traj) {
//   std::string path ="/home/ljin/application/Public/Muon-GM2-cc/jobs/convert/data/24D/light-minus-heavy-1/traj=0" + std::to_string(traj); 
//   assert(dirExists(path));
//   return path;
// }


std::string three_point_24ID(int traj) {
  std::string path = "/gpfs/mira-fs1/projects/HadronicLight_4/ctu/hlbl/hlbl-pion/TwoPointWallCorrField/24D-0.00107/ama/t-min=0010/results=" + std::to_string(traj) + "/avg ; type=0";    
  assert(dirExists(path));
  return path;
}

std::string three_point_disc_24ID(int traj) {
  std::string path = "/projects/CSC249ADSE03/yidizhao/pGG_config/24ID/disc/t-min=20/pGG_disc." + std::to_string(traj);
  // std::string path = "/projects/CSC249ADSE03/yidizhao/pGG_config/24ID/disc/t-min=10/pGG_disc." + std::to_string(traj);
  assert(dirExists(path));
  return path;
}

std::string three_point_disc2_32ID(int traj) {
  std::string path = "/gpfs/mira-fs0/projects/CSC249ADSE03/yidizhao/pGG_config/32ID/disc_2/pGG_disc." + std::to_string(traj);
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


std::string three_point_path_32IDF(int traj) {

  std::string path = "/home/ljin/application/Public/Qlat-CPS-cc/jobs/em-corr/results/32Dfine-0.0001/results="+ std::to_string(traj) + "/contraction-with-point/pion_gg/decay_cheng";
  assert(dirExists(path));
	return path;
}

std::string three_point_48I(int traj) {

  std::string path = "/home/ljin/application/Public/Qlat-CPS-cc/jobs/em-corr/results/48I-0.00078/results=" + std::to_string(traj) + "/contraction-with-point/pion_gg/decay_cheng";
  assert(dirExists(path));
	return path;
}


