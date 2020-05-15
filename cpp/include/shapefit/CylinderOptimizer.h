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
#include "objects/Cylinder.h"
#include "LocalParamGraff13.h"

/**
 * residual for Ceres optimization of cylinder parameters
 */
struct CylinderResidual {

    CylinderResidual(Eigen::Vector3f p) :
        p_(p.cast<double>())
    {}

    template <typename T>
    bool operator()(const T* ca, const T* r, T* residual) const {
        
        const T dx = ca[0] - p_[0];
        const T dy = ca[1] - p_[1];
        const T dz = ca[2] - p_[2];
        const T dTa = dx * ca[3] + dy * ca[4] + dz * ca[5];
        
        const T dist_to_center = ceres::sqrt(dx * dx + dy * dy + dz * dz - dTa * dTa);
        
        residual[0] = dist_to_center - r[0];
        
        return true;    
    }

    private:
        Eigen::Vector3d p_;
};

/*
 * actual Cylinder optimizer class
 */
template <class Pcd>
class CylinderOptimizer : public ShapeOptimizer<Pcd> {

    Cylinder<double>* C_;
    
    virtual void add_residual(Eigen::Vector3f point) {   
        CylinderResidual* res = new CylinderResidual(point);
        ceres::CostFunction* cost = new ceres::AutoDiffCostFunction<CylinderResidual, 1, 6, 1>(res);
        this->problem_->AddResidualBlock(cost, new TruncatedLoss(this->scale_), C_->data(), &(C_->data()[6]));  
    }
    
    virtual void set_param() {
        this->problem_->SetParameterization(C_->data(), new LocalParamGraff13);
    }
    
public:

    CylinderOptimizer(Cylinder<double>* C) :
        C_(C)
    {}

};
