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
#include "detector/detectors.h"
#include "detection_metrics.h"

int detect(int argc, char** argv) {

    Timer T;

    std::string folder, infile, results, otype;
    int det_settings = 64;
    size_t num_points = 4096;
    bool is_scene = true, is_object = false;


    CLI::App app{"Detect objects in .ply point cloud"};
    app.add_option("--folder", folder, "folder of input point cloud")->required();
    app.add_option("--in", infile, "input point cloud, in .ply format")->required();
    app.add_option("--results", results, "folder where results shall be written")->required();
    app.add_option("--type", otype, "type of object to be detected")->required();
    app.add_option("--settings", det_settings, "int specifying DetectorSettings");
    app.add_option("--np", num_points, "number of input points to detection stage");
    app.add_flag("--o", is_object, "flag for object parts detection");

    try {
        app.parse(argc, argv);
    } catch (const CLI::ParseError& e) {
        return app.exit(e);
    }

    if (is_object) {
        is_scene = false;
    }

    std::vector<int> options_vec;
    if (det_settings < 64) {
        options_vec.push_back(det_settings);
    }
    else {
        options_vec = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 20, 21, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 52, 53};
    }

    enum class Object {
        ALL,
        PLANE,
        SPHERE,
        CYLINDER,
        CONE,
        PLANE_SPHERE,
        PLANE_CYL,
        PLANE_CONE
    };
    Object O;

    if (otype == "plane") {
        O = Object::PLANE;
    }
    else if (otype == "sphere") {
        O = Object::SPHERE;
    }
    else if (otype == "cylinder") {
        O = Object::CYLINDER;
    }
    else if (otype == "cone") {
        O = Object::CONE;
    }
    else if (otype == "plane-sphere") {
        O = Object::PLANE_SPHERE;
    }
    else if (otype == "plane-cyl") {
        O = Object::PLANE_CYL;
    }
    else if (otype == "plane-cone") {
        O = Object::PLANE_CONE;
    }
    else if (otype == "all") {
        O = Object::ALL;
    }
    else {
        std::cerr << "Your specified object type is not supported (yet)." << std::endl;
        return 1;
    }

    T.tic();
    PointCloudData pc(folder + infile);
    T.toc("Load ply file");

    std::cout << "Number of points: " << pc.num_points() << std::endl;
    std::cout << "Diameter of scene: " << pc.diameter() << std::endl;

    T.tic();
    PcSampler<PointCloudData, float> pcs;
    std::vector<size_t> idx = pcs.sample_random(pc, num_points);
    pc.sample(idx);
    T.toc("Downsample point cloud");

    std::cout << "Number of points: " << pc.num_points() << std::endl;
    std::cout << "Diameter of scene: " << pc.diameter() << std::endl;

    PpfDetector<PointCloudData, float, float>* detector;

    for (const int& o : options_vec) {
        std::bitset<6> options(o);
        DetectorSettings settings(options);
        std::cout << std::endl << ".........." << std::endl;
        settings.print();
        std::cout << std::endl;
        switch (O) {
            case Object::PLANE :
                detector = new PlaneDetector<PointCloudData, float, float>(pc.diameter(), is_scene, settings);
                break;
            case Object::SPHERE :
                detector = new SphereDetector<PointCloudData, float, float>(pc.diameter(), is_scene, settings);
                break;
            case Object::CYLINDER :
                detector = new CylinderDetector<PointCloudData, float, float>(pc.diameter(), is_scene, settings);
                break;
            case Object::CONE :
                detector = new ConeDetector<PointCloudData, float, float>(pc.diameter(), is_scene, settings);
                break;
            case Object::PLANE_SPHERE :
                detector = new PlaneSphereDetector<PointCloudData, float, float>(pc.diameter(), is_scene, settings);
                break;
            case Object::PLANE_CYL :
                detector = new PlaneCylinderDetector<PointCloudData, float, float>(pc.diameter(), is_scene, settings);
                break;
            case Object::PLANE_CONE :
                detector = new PlaneConeDetector<PointCloudData, float, float>(pc.diameter(), is_scene, settings);
                break;
            case Object::ALL :
                detector = new PrimitiveDetector<PointCloudData, float, float>(pc.diameter(), is_scene, settings);
                break;
            default:
                std::cerr << "Your specified object type is not supported (yet)." << std::endl;
                return 1;
        }

        T.tic();
        detector->cast_votes(pc);
        T.toc("Voting");
        detector->print_info();
//        detector->print_candidates();
        T.tic();
        detector->cluster_candidates();
        detector->cluster_candidates();
        T.toc("Two rounds of clustering");
        detector->print_info();
        detector->print_candidates(1e3);

        detector->write_results(results, infile, 12, .35*num_points);

        std::cout << "# inliers, coverage, 80th percentiles, mean residuals and RMS values for candidate objects:" << std::endl;
        std::vector<MyObject<float>*> candidates = detector->candidate_vector(12, .35*num_points);
        for (auto& obj : candidates) {
            std::cout << metric_inliers<PointCloudData, float>(pc, obj) << "\t"
                      << metric_coverage<PointCloudData, float>(pc, obj)*100 << "\%\t"
                      << metric_p80<PointCloudData, float>(pc, obj) << "\t"
                      << metric_mean<PointCloudData, float>(pc, obj) << "\t"
                      << metric_rms<PointCloudData, float>(pc, obj) << std::endl;
        }
    }

    delete detector;

    return 0;
}
