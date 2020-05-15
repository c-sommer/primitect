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

#include "ShapeOptimizer.h"
#include "objects/Sphere.h"

/**
 * residual for Ceres optimization of sphere parameters
 */
struct SphereResidual {

    SphereResidual(Eigen::Vector3f p) :
        p_(p.cast<double>())
    {}

    template <typename T>
    bool operator()(const T* cr, T* residual) const {
        
        const T dx = cr[0] - p_[0];
        const T dy = cr[1] - p_[1];
        const T dz = cr[2] - p_[2];
        
        const T dist_to_center = ceres::sqrt(dx * dx + dy * dy + dz * dz);
        
        residual[0] = dist_to_center - cr[3];
        
        return true;    
    }

    private:
        Eigen::Vector3d p_;
};

/*
 * actual sphere optimizer class
 */
template <class Pcd>
class SphereOptimizer : public ShapeOptimizer<Pcd> {

    Sphere<double>* S_;
    
    virtual void add_residual(Eigen::Vector3f point) {   
        SphereResidual* res = new SphereResidual(point);
        ceres::CostFunction* cost = new ceres::AutoDiffCostFunction<SphereResidual, 1, 4>(res);
        this->problem_->AddResidualBlock(cost, new TruncatedLoss(this->scale_), S_->data());    
    }
    
public:

    SphereOptimizer(Sphere<double>* S) :
        S_(S)
    {}

};
