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
#include "objects/Plane.h"
#include "LocalParamS2.h"

/**
 * residual for Ceres optimization of Plane parameters
 */
struct PlaneResidual {

    PlaneResidual(Eigen::Vector3f p) :
        p_(p.cast<double>())
    {}

    template <typename T>
    bool operator()(const T* n, const T* d, T* residual) const {

        residual[0] = n[0] * p_[0] + n[1] * p_[1] + n[2] * p_[2] + d[0];
        
        return true;    
    }

    private:
        Eigen::Vector3d p_;
};

/*
 * actual Plane optimizer class
 */
template <class Pcd>
class PlaneOptimizer : public ShapeOptimizer<Pcd> {

    Plane<double>* P_;
    
    virtual void add_residual(Eigen::Vector3f point) {   
        PlaneResidual* res = new PlaneResidual(point);
        ceres::CostFunction* cost = new ceres::AutoDiffCostFunction<PlaneResidual, 1, 3, 1>(res);
        this->problem_->AddResidualBlock(cost, new TruncatedLoss(this->scale_), P_->data(), &(P_->data()[3]));  
    }
    
    virtual void set_param() {
        this->problem_->SetParameterization(P_->data(), new LocalParamS2);
    }
    
public:

    PlaneOptimizer(Plane<double>* P) :
        P_(P)
    {}

};
