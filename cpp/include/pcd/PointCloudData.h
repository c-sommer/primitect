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

#ifndef POINT_CLOUD_DATA_H_
#define POINT_CLOUD_DATA_H_

#include <vector>
#include <iostream>
#include <Eigen/Dense>

class PointCloudData {

    std::vector<Eigen::Vector3f> points_;
    std::vector<Eigen::Vector3f> normals_;
    size_t num_points_;

public:

    PointCloudData(std::vector<Eigen::Vector3f>& points, std::vector<Eigen::Vector3f>& normals) :
        points_(points),
        normals_(normals),
        num_points_(points.size())
    {
    }

    PointCloudData(std::string filepath);

    float point(size_t idx, size_t dim) const {
        return points_[idx][dim];
    }

    float normal(size_t idx, size_t dim) const {
        return normals_[idx][dim];
    }

    Eigen::Vector3f point(size_t idx) const {
        return points_[idx];
    }

    Eigen::Vector3f normal(size_t idx) const {
        return normals_[idx];
    }

    size_t num_points() const {
        return num_points_;
    }

    void normalize() {
        for (auto& n : normals_) {
            n.normalize();
        }
    }

    void sample(std::vector<size_t> indices) {
        std::vector<Eigen::Vector3f> new_points, new_normals;
        for (const auto& idx : indices) {
            new_points.push_back(points_[idx]);
            new_normals.push_back(normals_[idx]);
        }
        points_ = new_points;
        normals_ = new_normals;
        num_points_ = points_.size();
    }

    void scale(float scale) {
        for (auto& p : points_) {
            p *= scale;
        }
    }

    float diameter() const {
        float min_x = std::numeric_limits<float>::max(),    min_y = min_x, min_z = min_y;
        float max_x = std::numeric_limits<float>::lowest(), max_y = max_x, max_z = max_y;
        for (auto& p : points_) {
            if (p[0] < min_x) min_x = p[0];
            if (p[0] > max_x) max_x = p[0];
            if (p[1] < min_y) min_y = p[1];
            if (p[1] > max_y) max_y = p[1];
            if (p[2] < min_z) min_z = p[2];
            if (p[2] > max_z) max_z = p[2];
        }
        const float dx = max_x - min_x;
        const float dy = max_y - min_y;
        const float dz = max_z - min_z;
        return std::sqrt(dx * dx + dy * dy + dz * dz);
    }

    void write(std::string filepath);

};

#endif // POINT_CLOUD_DATA_H_
