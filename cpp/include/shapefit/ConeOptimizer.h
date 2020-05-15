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
#include "objects/Cone.h"
#include "LocalParamS2.h"
#include "LocalParamTheta.h"

/**
 * cost function for Ceres optimization of Cone parameters
 * due to if-condition in residual computation, automatic differentiation is not straightforward
 */
class ConeCostFunction : public ceres::SizedCostFunction<1, 3, 3, 3> {

    Eigen::Vector3d p_;

public:

    ConeCostFunction(Eigen::Vector3f p) :
        ceres::SizedCostFunction<1, 3, 3, 3>(),
        p_(p.cast<double>())
    {}

    virtual ~ConeCostFunction() {}
    
    virtual bool Evaluate(double const* const* parameters,
                          double* residual,
                          double** jacobians) const {
                          
        const double* c = parameters[0];
        const double* a = parameters[1];
        const double* theta3 = parameters[2];
        
        const double dx = p_[0] - c[0];
        const double dy = p_[1] - c[1];
        const double dz = p_[2] - c[2];
        
        const double d_sq = dx * dx + dy * dy + dz * dz;
        
        const double h = dx * a[0] + dy * a[1] + dz * a[2];
        const double r = std::sqrt(d_sq - h * h);
        
        if (h * theta3[1] + r * theta3[0] < 0.) {
            residual[0] = - std::sqrt(d_sq);
            if (jacobians && jacobians[0]) {
                // derivatives w.r.t. c
                const double d_inv = - 1. / residual[0];
                jacobians[0][0] = dx * d_inv;
                jacobians[0][1] = dy * d_inv;
                jacobians[0][2] = dz * d_inv;
                // derivatives w.r.t. a are zero
                jacobians[1][0] = jacobians[1][1] = jacobians[1][2] = 0.;
                // derivatives w.r.t. theta3 are zero
                jacobians[2][0] = jacobians[2][1] = jacobians[2][2] = 0.;
            }
        }
        else {
            residual[0] = h * theta3[0] - r * theta3[1];
            if (jacobians && jacobians[0]) {
                const double r_inv = 1. / r;
                const double coeff = theta3[0] + h * r_inv * theta3[1];
                // derivatives w.r.t. c
                jacobians[0][0] = - coeff * a[0] + theta3[1] * r_inv * dx;
                jacobians[0][1] = - coeff * a[1] + theta3[1] * r_inv * dy;
                jacobians[0][2] = - coeff * a[2] + theta3[1] * r_inv * dz;
                // derivatives w.r.t. a
                jacobians[1][0] = coeff * dx;
                jacobians[1][1] = coeff * dy;
                jacobians[1][2] = coeff * dz;
                // derivatives w.r.t. theta3
                jacobians[2][0] = h;
                jacobians[2][1] = -r;
                jacobians[2][2] = 0.;
            }
        }

        return true;
    }
};




/*
 * actual Cone optimizer class
 */
template <class Pcd>
class ConeOptimizer : public ShapeOptimizer<Pcd> {

    Cone<double>* C_;
    
    virtual void add_residual(Eigen::Vector3f point) {   
//        ConeResidual* res = new ConeResidual(point);
//        ceres::CostFunction* cost = new ceres::AutoDiffCostFunction<ConeResidual, 1, 3, 3, 3>(res);
        ceres::CostFunction* cost = new ConeCostFunction(point);
        double* params = C_->data();
        this->problem_->AddResidualBlock(cost, new TruncatedLoss(this->scale_), &(params[0]), &(params[3]), &(params[6]));
    }
    
    virtual void set_param() {
        this->problem_->SetParameterization(&(C_->data()[3]), new LocalParamS2);
        this->problem_->SetParameterization(&(C_->data()[6]), new LocalParamTheta);
    }
    
public:

    ConeOptimizer(Cone<double>* C) :
        C_(C)
    {}

};
