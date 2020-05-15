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

#include "PointCloudData.h"
#include <fstream>
#include <tinyply/tinyply.h>

PointCloudData::PointCloudData(std::string filepath) {

	try {
	
		std::ifstream ss(filepath, std::ios::binary);
		if (ss.fail())
		    throw std::runtime_error("failed to open " + filepath);

		tinyply::PlyFile file;
		file.parse_header(ss);

		// Tinyply treats parsed data as untyped byte buffers. See below for examples.
		std::shared_ptr<tinyply::PlyData> points, normals;

		try {
		    points = file.request_properties_from_element("vertex", { "x", "y", "z" });
		} catch (const std::exception & e) {
		    std::cerr << "tinyply exception: " << e.what() << std::endl;
		}

		try {
		    normals = file.request_properties_from_element("vertex", { "nx", "ny", "nz" });
		} catch (const std::exception & e) {
		    std::cerr << "tinyply exception: " << e.what() << std::endl;
		}

		file.read(ss);

		// type casting to your own native types
		points_ = std::vector<Eigen::Vector3f>(points->count);
		std::memcpy(points_.data(), points->buffer.get(), points->buffer.size_bytes());
		normals_ = std::vector<Eigen::Vector3f>(normals->count);
		std::memcpy(normals_.data(), normals->buffer.get(), normals->buffer.size_bytes());
		
		num_points_ = points_.size();
		
	} catch (const std::exception & e) {
	
		std::cerr << "Caught tinyply exception: " << e.what() << std::endl;
	}
	normalize();
}

void PointCloudData::write(std::string filepath) {

    std::filebuf fb_binary;
	fb_binary.open(filepath, std::ios::out | std::ios::binary);
    std::ostream outstream_binary(&fb_binary);
	if (outstream_binary.fail())
	    throw std::runtime_error("failed to open " + filepath);

    tinyply::PlyFile ply_file;
    ply_file.add_properties_to_element("vertex", { "x", "y", "z" }, tinyply::Type::FLOAT32, num_points_, reinterpret_cast<uint8_t*>(points_.data()), tinyply::Type::INVALID, 0);
    ply_file.add_properties_to_element("vertex", { "nx", "ny", "nz" }, tinyply::Type::FLOAT32, num_points_, reinterpret_cast<uint8_t*>(normals_.data()), tinyply::Type::INVALID, 0);
	ply_file.write(outstream_binary, true);

}
