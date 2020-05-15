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

#include "pcd/PointCloudData.h"
#include "objects/MyObject.h"

template <class Pcd, typename T>
size_t metric_inliers(const Pcd& pc, MyObject<T>* O, const T threshold = 0.05) {

    size_t counter = 0;
    size_t N = pc.num_points();

    for (size_t idx=0; idx<N; ++idx) {
        T psi = O->sdf(pc.point(idx).template cast<T>());
        if (std::abs(psi) < threshold) {
            T weight = pc.normal(idx).template cast<T>().dot(O->normal(pc.point(idx).template cast<T>()));
            if (weight * weight > .95) {
                ++counter;
            }
        }
    }
    
    return counter;

}

template <class Pcd, typename T>
T metric_coverage(const Pcd& pc, MyObject<T>* O, const T threshold = 0.05) {
    
    return metric_inliers(pc, O, threshold) / static_cast<T>(pc.num_points());

}

template <class Pcd, typename T>
T metric_rms(const Pcd& pc, MyObject<T>* O, const T threshold = 0.05) {

    T sum = 0.;
    size_t counter = 0;
    size_t N = pc.num_points();

    for (size_t idx=0; idx<N; ++idx) {
        T psi = O->sdf(pc.point(idx).template cast<T>());
        if (psi > -threshold && psi < threshold) {
            T weight = pc.normal(idx).template cast<T>().dot(O->normal(pc.point(idx).template cast<T>()));
            if (weight * weight > .95) {
                sum += psi * psi;
                ++counter;
            }
        }
    }
    
    if (counter<1) {
        return 0;
    }
    
    return std::sqrt(sum / counter);

}

template <class Pcd, typename T>
T metric_mean(const Pcd& pc, MyObject<T>* O, const T threshold = 0.05) {

    T sum = 0.;
    size_t counter = 0;
    size_t N = pc.num_points();

    for (size_t idx=0; idx<N; ++idx) {
        T psi = std::abs(O->sdf(pc.point(idx).template cast<T>()));
        if (psi < threshold) {
            T weight = pc.normal(idx).template cast<T>().dot(O->normal(pc.point(idx).template cast<T>()));
            if (weight * weight > .95) {
                sum += psi;
                ++counter;
            }
        }
    }
    
    if (counter<1) {
        return 0.;
    }
    
    return sum / counter;

}

template <class Pcd, typename T>
T metric_p80(const Pcd& pc, MyObject<T>* O, const T threshold = 0.05) {

    T sum = 0.;
    std::vector<T> residuals;
    size_t N = pc.num_points();

    for (size_t idx=0; idx<N; ++idx) {
        T psi = std::abs(O->sdf(pc.point(idx).template cast<T>()));
        if (psi < threshold) {
            T weight = pc.normal(idx).template cast<T>().dot(O->normal(pc.point(idx).template cast<T>()));
            if (weight * weight > .95) {
                residuals.push_back(psi);
            }
        }
    }
    
    if (residuals.size()<1) {
        return 0.;
    }
    
    size_t idx = static_cast<size_t>(.80 * residuals.size());
    std::nth_element(residuals.begin(), residuals.begin() + idx, residuals.end());
    
    return residuals[idx];
}

template <class Pcd, typename T>
void print_metrics(const Pcd& pc, MyObject<T>* O, const T threshold = 0.05) {

    std::cout << "Inliers:\t" << metric_inliers<Pcd, T>(pc, O) << "\t"
              << "Coverage:\t" << metric_coverage<Pcd, T>(pc, O) * 100 << "\%\t"
              << "80th percentile:\t" << metric_p80<Pcd, T>(pc, O) << "\t"
              << "mean:\t" << metric_mean<Pcd, T>(pc, O) << "\t"
              << "RMS:\t" << metric_rms<Pcd, T>(pc, O) << std::endl;

    return;
}
