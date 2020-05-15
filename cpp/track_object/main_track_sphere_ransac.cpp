/**
GPL (v3+) License

This file is part of the code accompanying the paper
PrimiTect: Fast Continuous Hough Voting for Primitive Detection
by C. Sommer, Y. Sun, E. Bylow and D. Cremers,
accepted for publication in the IEEE International Conference on Robotics and Automation (ICRA) 2020.

Copyright (c) 2019, Christiane Sommer.
Adapted and extended from CGAL Shape Detection examples
  (c) 2015, Sven Oesau, Yannick Verdie, Clément Jamin, Pierre Alliez.
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

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/IO/read_xyz_points.h>
#include <CGAL/Point_with_normal_3.h>
#include <CGAL/Shape_detection_3.h>
#include <CGAL/Timer.h>
#include <CGAL/number_utils.h>
#include <CGAL/property_map.h>
#include <CGAL/random_simplify_point_set.h>
#include <CLI/CLI.hpp>
#include <fstream>
#include <iostream>


// Type declarations
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::FT FT;
typedef std::pair<Kernel::Point_3, Kernel::Vector_3> Point_with_normal;
typedef std::vector<Point_with_normal> Pwn_vector;
typedef CGAL::First_of_pair_property_map<Point_with_normal> Point_map;
typedef CGAL::Second_of_pair_property_map<Point_with_normal> Normal_map;

// In Shape_detection_traits the basic types, i.e., Point and Vector types
// as well as iterator type and property maps, are defined.
typedef CGAL::Shape_detection_3::Shape_detection_traits<Kernel,Pwn_vector, Point_map, Normal_map> Traits;
typedef CGAL::Shape_detection_3::Efficient_RANSAC<Traits> Efficient_ransac;
typedef CGAL::Shape_detection_3::Plane<Traits> Plane;
typedef CGAL::Shape_detection_3::Sphere<Traits> Sphere;
typedef CGAL::Shape_detection_3::Cylinder<Traits> Cylinder;
typedef CGAL::Shape_detection_3::Cone<Traits> Cone;

int SPHERE_TYPE = 1;
int PLANE_TYPE = 2;
int CYLINDER_TYPE = 3;
int CONE_TYPE = 4;


void label_type(Eigen::VectorXi& label_per_point, Eigen::VectorXi& type_per_point, int instance, int type, const std::vector<size_t>& indices) {
  for (int i = 0 ; i < indices.size(); i++){
    label_per_point(indices[i]) = instance;
    type_per_point(indices[i]) = type;
  }
}

int main(int argc, char** argv) {

    std::string folder, results, pc_file;
    std::string otype = "sphere";
    size_t num_points = 2048;
    
    CLI::App app{"Track sphere in .ply point clouds"};
    app.add_option("--folder", folder, "folder of input point cloud")->required();
    app.add_option("--results", results, "folder where results shall be written")->required();
    app.add_option("--np", num_points, "number of input points to detection stage");

    try {
        app.parse(argc, argv);
    } catch (const CLI::ParseError& e) {
        return app.exit(e);
    }   

  // Points with normals.
  Pwn_vector points;

  // Sets parameters for shape detection.
  Efficient_ransac::Parameters parameters;
  // Sets probability to miss the largest primitive at each iteration. (Default: 5%)
  // parameters.probability = 0.05;
  // Set minimium number of points that a shape must contain. (Default: 1% of points)
  // parameters.min_points = 200;
  // Sets maximum Euclidean distance between a point and a shape. (Default: 1% of BB diagonal)
  parameters.epsilon = 0.05;
  // Sets maximum Euclidean distance between points to be clustered. (Default: 1% of BB diagonal)
  // parameters.cluster_epsilon = 0.01;
  // Sets maximum normal deviation as cos of surface normal and point normal. (Default: .9 = cos(25.8°))
  // 0.9 < dot(surface_normal, point_normal);
  parameters.normal_threshold = std::cos(10. * M_PI / 180.);
  
  
    std::ifstream infile(folder + "files_xyz.txt");
    std::ofstream outfile(results + "spheres.txt");
    
    while (!infile.eof()) {
    
        infile >> pc_file;

        // Loads point set from a file.
        // read_xyz_points_and_normals takes an OutputIterator for storing the points
        // and a property map to store the normal vector with each point.
        std::ifstream stream(folder + pc_file);
        if (!stream || !CGAL::read_xyz_points(stream, std::back_inserter(points), CGAL::parameters::point_map(Point_map()).normal_map(Normal_map()))) {
            std::cerr << "Error: cannot read file" << std::endl;
            return EXIT_FAILURE;
        }

        if (num_points > 0) {
            points.erase(CGAL::random_simplify_point_set(points, 100.*(1. - double(num_points) / double(points.size()))), points.end());
        }

        // Instantiates shape detection engine.
        Efficient_ransac ransac;
        // Provides the input data.
        ransac.set_input(points);
        // Registers detection of objects
        ransac.add_shape_factory<Sphere>();
        // Build internal data structures.
        ransac.preprocess();
        // Perform detection
        Efficient_ransac::Shape_range shapes = ransac.shapes();
        ransac.detect(parameters);
        shapes = ransac.shapes();

        // Write first detection to file
        Efficient_ransac::Shape_range::iterator it = shapes.begin();
        Sphere* sph = dynamic_cast<Sphere*>(it->get());
        outfile << sph->center() << " " << sph->radius() << std::endl;      

    }
    
    infile.close();
    outfile.close();

  return EXIT_SUCCESS;
}
