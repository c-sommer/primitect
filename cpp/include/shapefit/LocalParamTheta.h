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

#ifndef LOCAL_PARAM_THETA_H_
#define LOCAL_PARAM_THETA_H_

#include <ceres/local_parameterization.h>

class LocalParamTheta : public ceres::LocalParameterization {

    static constexpr double PI_2 = .5 * M_PI;

public:

    LocalParamTheta() {}
    
    virtual ~LocalParamTheta() {}

    // Theta plus operation for Ceres
    //
    virtual bool Plus(double const* theta3, double const* delta,
                      double* theta3_plus_delta) const {
                      
        const double theta_plus_delta = theta3[2] + delta[0];
        // TODO: possibly truncate theta to [0, pi/2]
        
        theta3_plus_delta[0] = std::sin(theta_plus_delta);
        theta3_plus_delta[1] = std::sqrt(1. - theta3_plus_delta[0] * theta3_plus_delta[0]);
        theta3_plus_delta[2] = theta_plus_delta;

        return true;

    }

    // Jacobian of Theta plus operation for Ceres
    //
    virtual bool ComputeJacobian(double const* theta3,
                                 double* J) const {
        
        // Jacobian has size 3x1
        J[0] = theta3[1];
        J[1] = -theta3[0];
        J[2] = 1.;
        
        return true;

    }

    virtual int GlobalSize() const { return 3; }

    virtual int LocalSize() const { return 1; }
};

#endif // LOCAL_PARAM_THETA_H_
