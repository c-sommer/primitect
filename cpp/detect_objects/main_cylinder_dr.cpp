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
#include "detector/CylinderDetector.h"
#include "shapefit/CylinderOptimizer.h"
#include "detection_metrics.h"

int main(int argc, char** argv) {
    
    Timer T;
    
    std::string folder, infile;

    CLI::App app{"Detect objects in .ply point cloud"};
    app.add_option("--folder", folder, "folder of input and output point cloud")->required();
    app.add_option("--in", infile, "input point cloud, in .ply format")->required();

    try {
        app.parse(argc, argv);
    } catch (const CLI::ParseError& e) {
        return app.exit(e);
    }
    
    T.tic();
    PointCloudData pc(folder + infile);
    T.toc("Load ply file");
    
    DetectorSettings settings(std::bitset<6>(39));
    CylinderDetector<PointCloudData, float, float>* detector = new CylinderDetector<PointCloudData, float, float>(settings);
    CylinderOptimizer<PointCloudData>* optimizer = nullptr;
    
    T.tic();
    detector->cast_votes(pc);
    T.toc("Voting");
    T.tic();
    detector->cluster_candidates();
    detector->cluster_candidates();
    T.toc("Two rounds of clustering");
    
    for (auto el : detector->candidates()) {
    
        if (el.first < 1e3) break;

        Cylinder<double> C = el.second->cast<double>();
        MyObject<double>* obj = static_cast<MyObject<double>*>(&C);
                          
        std::cout << std::endl << ".......... # votes: " << el.first << std::endl;
        std::cout << "Before refinement:" << std::endl;
        std::cout << C << std::endl;
        print_metrics<PointCloudData, double>(pc, obj);
    
        optimizer = new CylinderOptimizer<PointCloudData>(&C);
        T.tic();
        bool conv = optimizer->optimize(pc);
        T.toc("Parameter refinement");
        if (conv) {
            obj = static_cast<MyObject<double>*>(&C);
            std::cout << "After refinement:" << std::endl;
            std::cout << C << std::endl;
            print_metrics<PointCloudData, double>(pc, obj);
        }
        else {
            std::cout << "No convergence." << std::endl;
        }

    }
        
    delete detector;
    if (optimizer) delete optimizer;
    
    return 0;
}
