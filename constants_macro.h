#pragma once

#define TIME_LIMIT 16
#define SPACE_LIMIT 16

// const std::string SLOPPY_WALL_PATH = "";

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


// // sloppy
//
// // const std::string point_file = "/home/ljin/application/Public/Muon-GM2-cc/jobs/" + ENSEMBLE + ssprintf("/discon-1/results/results=%d/checkpoint/computeContractionInf", i);
// // const std::string wall_file = "/home/ljin/application/Public/Qlat-CPS-cc/jobs/" + ENSEMBLE + ssprintf("/wall-src/results/32D-0.00107/results=%d/checkpoint.txt", i);
//
// const std::string WALL_SRC_PATH = "/home/ljin/application/Public/Qlat-CPS-cc/jobs/" + ENSEMBLE + ssprintf("/wall-src/results/32D-0.00107/results=%d/huge-data/wall_src_propagator", traj_list[i]);
// const std::string GAUGE_TRANSFORM_PATH = "/home/ljin/application/Public/Qlat-CPS-cc/jobs/" + ENSEMBLE + ssprintf("/wall-src/results/32D-0.00107/results=%d/huge-data/gauge-transform", traj_list[i]);
// const std::string POINT_SRC_PATH = "/home/ljin/application/Public/Muon-GM2-cc/jobs/" + ENSEMBLE + ssprintf("/discon-1/results/prop-hvp ; results=%d/huge-data/prop-point-src", traj_list[i]);
//
// // exact
//     // std::string point_file = "/home/ljin/application/Public/Muon-GM2-cc/jobs/" + ENSEMBLE + ssprintf("/discon-1/results/results=%d/checkpoint/computeContractionInf", i);
//     // std::string wall_file = "/home/ljin/application/Public/Qlat-CPS-cc/jobs/" + ENSEMBLE + ssprintf("/wall-src/results/32D-0.00107/results=%d/checkpoint.txt", i);
//     // std::string exact_wall_file = "/home/ljin/application/Public/Qlat-CPS-cc/jobs/" + ENSEMBLE + ssprintf("/wall-src-exact-2/results/32D-0.00107/results=%d/checkpoint.txt", i);
// 		//
// const std::string WALL_SRC_PATH = "/home/ljin/application/Public/Qlat-CPS-cc/jobs/" + ENSEMBLE + ssprintf("/wall-src/results/32D-0.00107/results=%d/huge-data/wall_src_propagator", traj_list[i]);
// const std::string EXACT_WALL_SRC_PATH = "/home/ljin/application/Public/Qlat-CPS-cc/jobs/" + ENSEMBLE + ssprintf("/wall-src-exact-2/results/32D-0.00107/results=%d/huge-data/wall_src_propagator", traj_list[i]);
// const std::string GAUGE_TRANSFORM_PATH = "/home/ljin/application/Public/Qlat-CPS-cc/jobs/" + ENSEMBLE + ssprintf("/wall-src/results/32D-0.00107/results=%d/huge-data/gauge-transform", traj_list[i]);
// const std::string POINT_SRC_PATH = "/home/ljin/application/Public/Muon-GM2-cc/jobs/" + ENSEMBLE + ssprintf("/discon-1/results/prop-hvp ; results=%d/huge-data/prop-point-src", traj_list[i]);
//
// // wall src for 24D wall-to-wall
//     const std::string WALL_SRC_PATH = ssprintf("/home/ljin/application/Public/Qlat-CPS-cc/jobs/24D/wall-src/results/results=%d/huge-data/wall_src_propagator", traj_list[i]);
//
// // wall src for 32D wall-to-wall This is for sloppy
//     std::string wall_file = ssprintf("/home/ljin/application/Public/Qlat-CPS-cc/jobs/32D/wall-src/results/32D-0.00107/results=%d/checkpoint.txt", traj);
//
//     const std::string WALL_SRC_PATH = ssprintf("/home/ljin/application/Public/Qlat-CPS-cc/jobs/32D/wall-src/results/32D-0.00107/results=%d/huge-data/wall_src_propagator", traj_list[i]);
