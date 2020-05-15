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

#ifndef PC_SAMPLER_H_
#define PC_SAMPLER_H_

#include <iostream>
#include <numeric>
#include <algorithm>
#include <random>
#include <Eigen/Dense>
// #include "Grid.h"

template <class Pcd, typename T>
class PcSampler {

    using Vec3 = Eigen::Matrix<T, 3, 1>;

    T stepsize_ = 1.; // sample every point
    
    T sampling_ratio_ = 0.01; // 1% of total diameter
    
    T grid_length_ = 0.01;

    std::vector<size_t> sample_uniform(Pcd& pc);

    std::vector<size_t> sample_grid(Pcd& pc);
    
    void bounding_box(const Pcd& pc, Vec3& min_range, Vec3& max_range);
    
    std::vector<size_t> sort_indices(const Pcd& pc);

public:
    
    std::vector<size_t> sample_uniform(Pcd& pc, size_t num_sampled);
    
    std::vector<size_t> sample_uniform(Pcd& pc, float ratio) {
        set_stepsize(static_cast<T>(ratio));
        return sample_uniform(pc);
    }
    
    std::vector<size_t> sample_uniform(Pcd& pc, double ratio) {
        set_stepsize(static_cast<T>(ratio));
        return sample_uniform(pc);
    }
    
    std::vector<size_t> sample_random(Pcd& pc, size_t num_points);

    std::vector<size_t> sample_random(Pcd& pc, float ratio) {
        size_t num_points = ratio * pc.num_points();
        return sample_random(pc, num_points);
    }
    
    std::vector<size_t> sample_random(Pcd& pc, double ratio) {
        size_t num_points = ratio * pc.num_points();
        return sample_random(pc, num_points);
    }
    
    std::vector<size_t> sample_grid(Pcd& pc, T grid_length) {
        grid_length_ = grid_length;
        return sample_grid(pc);
    }
    
    std::vector<size_t> cutoff_z(Pcd& pc, size_t num_points);
    
    void set_stepsize(T ratio) {
        stepsize_ = T(1.) / ratio;
    }

};

template <class Pcd, typename T>
std::vector<size_t> PcSampler<Pcd, T>::sample_uniform(Pcd& pc) {
    T num_points_f = static_cast<T>(pc.num_points());
    stepsize_ = std::max(T(1.), std::min(num_points_f, stepsize_));
    
    std::vector<size_t> indices;
    for (T idx = 0; idx < pc.num_points(); idx += stepsize_) {
        indices.push_back(static_cast<size_t>(idx));
    }
    stepsize_ = 1.;
    return indices;
}

template <class Pcd, typename T>
std::vector<size_t> PcSampler<Pcd, T>::sample_uniform(Pcd& pc, size_t num_sampled) {
    T num_points_f = static_cast<T>(pc.num_points());
    stepsize_ = num_points_f / num_sampled;
    return sample_uniform(pc);
}

template <class Pcd, typename T>
std::vector<size_t> PcSampler<Pcd, T>::sample_random(Pcd& pc, size_t num_points) {

    std::vector<size_t> indices(pc.num_points());
    std::iota(indices.begin(), indices.end(), 0);
    
    if (pc.num_points() > num_points) {
        std::random_device rd;
        std::mt19937 g(rd());
        std::shuffle(indices.begin(), indices.end(), g);
        indices.resize(num_points);   
    }    
    return indices;
}

// For license reasons, we can unfortunately not provide the PcSampler<Pcd, T>::sample_grid(Pcd& pc) method in this repository.

template <class Pcd, typename T>
std::vector<size_t> PcSampler<Pcd, T>::sort_indices(const Pcd& pc) {

    // initialize original index locations
    std::vector<size_t> indices(pc.num_points());
    std::iota(indices.begin(), indices.end(), 0);
    
    // sort indices based on comparing z-values in point cloud pc
    std::sort(indices.begin(), indices.end(), [&pc](size_t i1, size_t i2) {return pc.point(i1, 2) < pc.point(i2, 2);});
    return indices;
}

template <class Pcd, typename T>
std::vector<size_t> PcSampler<Pcd, T>::cutoff_z(Pcd& pc, size_t num_points) {

    size_t N = pc.num_points();
    std::vector<size_t> indices(N);
    
    if (N <= num_points) {
        std::iota(indices.begin(), indices.end(), 0);
    }
    else {
        indices = sort_indices(pc);
        indices.resize(num_points);
    }
    
    return indices;
    
}

#endif // PC_SAMPLER_H_
