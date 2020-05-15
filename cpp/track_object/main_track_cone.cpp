/**
GPL (v3+) License

This file is part of the code accompanying the paper
PrimiTect: Fast Continuous Hough Voting for Primitive Detection
by C. Sommer, Y. Sun, E. Bylow and D. Cremers,
accepted for publication in the IEEE International Conference on Robotics and Automation (ICRA) 2020.

Copyright (c) 2019, Christiane Sommer.
All rights reserved.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

#include <iostream>
#include <bitset>
#include <Eigen/Dense>
#include <CLI/CLI.hpp>
#include "Timer.h"
#include "pcd/PointCloudData.h"
#include "sampling/PcSampler.hpp"
#include "detector/ConeDetector.h"
#include "detection_metrics.h"

int main(int argc, char** argv) {
    
    Timer T;
    
    std::string folder, results, pc_file;
    std::string otype = "cone";
    size_t num_points = 2048, s = 39;
    
    CLI::App app{"Track cone in .ply point clouds"};
    app.add_option("--folder", folder, "folder of input point cloud")->required();
    app.add_option("--results", results, "folder where results shall be written")->required();
    app.add_option("--np", num_points, "number of input points to detection stage");
    app.add_option("--settings", s, "settings");

    try {
        app.parse(argc, argv);
    } catch (const CLI::ParseError& e) {
        return app.exit(e);
    }   
    
    std::bitset<6> options(s);
    DetectorSettings settings(options);
    std::cout << std::endl << ".........." << std::endl;
    settings.print();
    std::cout << std::endl;
    
    PpfDetector<PointCloudData, float, float>* detector;
    
    std::ifstream infile(folder + "files.txt");
    std::ofstream outfile(results + "cones.txt");
    
    while (!infile.eof()) {
    
        infile >> pc_file;
    
        PointCloudData pc(folder + pc_file);
        PcSampler<PointCloudData, float> pcs;
        std::vector<size_t> idx = pcs.sample_random(pc, num_points);
        pc.sample(idx);
        
        detector = new ConeDetector<PointCloudData, float, float>(pc.diameter(), 1, settings);
        
        detector->cast_votes(pc);
        detector->cluster_candidates();
        detector->cluster_candidates();
        
        // TODO: only write first candidate to list!    
        std::vector<MyObject<float>*> candidates = detector->candidate_vector(12, .35*num_points);
        Cone<float>* cone = dynamic_cast<Cone<float>*>(candidates[0]);
        outfile << *cone << std::endl;
    }
    
    infile.close();
    outfile.close();
    delete detector;
    
    return 0;
}
