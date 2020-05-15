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

  // parse setttings from command line
  std::string folder, infile, results, otype;
  
  size_t num_points = 0;

  CLI::App app{"Detect objects in .xyz point cloud"};
  app.add_option("--folder", folder, "folder of input and output point cloud")->required();
  app.add_option("--in", infile, "input point cloud, in .xyz format")->required();
  app.add_option("--results", results, "folder where results shall be written")->required();
  app.add_option("--type", otype, "type of object to be detected")->required();
  app.add_option("--np", num_points, "number of input points to detection stage");

  try {
    app.parse(argc, argv);
  } catch (const CLI::ParseError& e) {
    return app.exit(e);
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

  // Points with normals.
  Pwn_vector points;

  // Measures time before loading the data.
  CGAL::Timer time;
  time.start();
  // Loads point set from a file.
  // read_xyz_points_and_normals takes an OutputIterator for storing the points
  // and a property map to store the normal vector with each point.
  std::ifstream stream(folder + infile);
  if (!stream || !CGAL::read_xyz_points(stream, std::back_inserter(points), CGAL::parameters::point_map(Point_map()).normal_map(Normal_map()))) {
    std::cerr << "Error: cannot read file" << std::endl;
    return EXIT_FAILURE;
  }
  // Measures time after preprocessing.
  time.stop();
  std::cout << "loading data took: " << time.time() * 1000 << "ms" << std::endl;
  
  if (num_points > 0) {
      points.erase(CGAL::random_simplify_point_set(points, 100.*(1. - double(num_points) / double(points.size()))), points.end());
  }

  // Instantiates shape detection engine.
  Efficient_ransac ransac;
  // Provides the input data.
  ransac.set_input(points);
  // Registers detection of objects
  switch (O) {
    case Object::PLANE :
      ransac.add_shape_factory<Plane>();
      break;
    case Object::SPHERE :
      ransac.add_shape_factory<Sphere>();
      break;
    case Object::CYLINDER :
      ransac.add_shape_factory<Cylinder>();
      break;
    case Object::CONE :
      ransac.add_shape_factory<Cone>();
      break;
    case Object::PLANE_SPHERE :
      ransac.add_shape_factory<Plane>();
      ransac.add_shape_factory<Sphere>();
      break;
    case Object::PLANE_CYL :
      ransac.add_shape_factory<Plane>();
      ransac.add_shape_factory<Cylinder>();
      break;
    case Object::PLANE_CONE :
      ransac.add_shape_factory<Plane>();
      ransac.add_shape_factory<Cone>();
      break;
    case Object::ALL :
      ransac.add_shape_factory<Plane>();
      ransac.add_shape_factory<Sphere>();
      ransac.add_shape_factory<Cylinder>();
      ransac.add_shape_factory<Cone>();
      break;
    default:
      std::cerr << "Your specified object type is not supported (yet)." << std::endl;
      return 1;
  }

  // Measures time before setting up the shape detection.
  time.start();
  // Build internal data structures.
  ransac.preprocess();
  // Measures time after preprocessing.
  time.stop();
  std::cout << "preprocessing took: " << time.time() * 1000 << "ms" << std::endl;

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

  // Perform detection several times and choose result with highest coverage.
  Efficient_ransac::Shape_range shapes = ransac.shapes();
  FT best_coverage = 0;
  for (size_t i = 0; i < 1; i++) {
    // Reset timer.
    time.reset();
    time.start();
    // Detects shapes.
    ransac.detect(parameters);
    // Measures time after detection.
    time.stop();
    // Compute coverage, i.e. ratio of the points assigned to a shape.
    FT coverage = FT(points.size() - ransac.number_of_unassigned_points())
                  / FT(points.size());
    // Prints number of assigned shapes and unsassigned points.
    std::cout << "time: " << time.time() * 1000 << "ms" << std::endl;
    std::cout << ransac.shapes().end() - ransac.shapes().begin() << " primitives, "
              << coverage << " coverage" << std::endl;

    // Choose result with highest coverage.
    if (coverage > best_coverage) {
      best_coverage = coverage;
      // Efficient_ransac::shapes() provides
      // an iterator range to the detected shapes.
      shapes = ransac.shapes();
    }
  }

  std::cout<<"number of instances: "<<shapes.size()<<std::endl;
  std::cout<<"point size N = "<<points.size()<<std::endl;
  Eigen::VectorXi label_per_point(points.size());
  label_per_point.setZero();

  Eigen::VectorXi type_per_point(points.size());
  type_per_point.setZero();

  std::string prefix = infile.substr(0, infile.size()-4);
  std::cout<<"file prefix: "<<prefix<<std::endl;

  std::ofstream labelsFile(results + prefix + "_points_with_labels.txt");
  if ( !labelsFile.is_open() ) {
    std::cout<<"label not open...\n";
    return -1;
  }
  std::ofstream planesFile(results + prefix + "_planes.txt");
  if ( !planesFile.is_open() ) {
    std::cout<<"plane not open...\n";
    return -1;
  }
  std::ofstream spheresFile(results + prefix + "_spheres.txt");
  if ( !spheresFile.is_open() ) {
    std::cout<<"sphere not open...\n";
    return -1;
  }
  std::ofstream cylindersFile(results + prefix + "_cylinders.txt");
  if ( !cylindersFile.is_open() ) {
    std::cout<<"cylinder not open...\n";
    return -1;

  }
  std::ofstream conesFile(results + prefix + "_cones.txt");
  if ( !conesFile.is_open() ) {
    std::cout<<"cone not open...\n";
    return -1;

  }

  int instance = 1;
  int type = 0;
  for (Efficient_ransac::Shape_range::iterator it = shapes.begin(); it != shapes.end(); it++) {
    
    if (Plane* plane = dynamic_cast<Plane*>(it->get()))
    {
      planesFile << instance << " " << plane->plane_normal() << " " << plane->d() << std::endl;
      type = 2;
    }
    else if (Cylinder* cyl = dynamic_cast<Cylinder*>(it->get()))
    {
      cylindersFile<<instance << " " << cyl->axis().point() << " " << cyl->axis().direction() << " " << cyl->radius() << std::endl;
      type = 3;
    }
    else if (Sphere* sph = dynamic_cast<Sphere*>(it->get()))
    {
      spheresFile << instance << " " << sph->center() << " " << sph->radius() << std::endl;
      type = 1;
    }
    else if (Cone* con = dynamic_cast<Cone*>(it->get()))
    {
      conesFile << instance << " " << con->apex() << " " << con->axis() << " " << con->angle() << std::endl;
      type = 4;
    }
    else
    {
      std::cout << (*it)->info() << std::endl;
    }
  
    // Iterates through point indices assigned to each detected shape.
    for (auto index_it = (*it)->indices_of_assigned_points().begin(); index_it != (*it)->indices_of_assigned_points().end(); index_it++) {      
      // Retrieves point
      const Point_with_normal &p = *(points.begin() + (*index_it));
      // write point and label to file
      labelsFile << p.first.x() << " " << p.first.y() << " " << p.first.z() << " " << instance << " " << type << std::endl;
    }

    instance++;
  }

  labelsFile.close();
  planesFile.close();
  spheresFile.close();
  cylindersFile.close();
  conesFile.close();

  return EXIT_SUCCESS;
}
