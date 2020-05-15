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

#ifndef LOCAL_PARAM_GRAFF13_H_
#define LOCAL_PARAM_GRAFF13_H_

#include <ceres/local_parameterization.h>

class LocalParamGraff13 : public ceres::LocalParameterization {

public:

    LocalParamGraff13() {}
    
    virtual ~LocalParamGraff13() {}

    // Graff13 plus operation for Ceres
    //
    virtual bool Plus(double const* ca, double const* delta,
                      double* ca_plus_delta) const { // ca: c (point on axis), a (axis direction)
        
        const double norm_inv = 1. / std::sqrt(1. + delta[2]*delta[2] + delta[3]*delta[3]);
        ca_plus_delta[0] =  ca[0] - ca[4] * delta[0] - ca[5] * delta[1];
        ca_plus_delta[3] = (ca[3] - ca[4] * delta[2] - ca[5] * delta[3]) * norm_inv;
        
        const double a_x_dc = (ca[5] * delta[0] - ca[4] * delta[1]);
        const double a_x_da = (ca[5] * delta[2] - ca[4] * delta[3]);
        
        if (ca[3] > 0.) { // ca[3] = a[x]
            const double ax1_inv = 1. / (1. + ca[3]);
            ca_plus_delta[1] =  ca[1] + ca[3] * delta[0] + ax1_inv * a_x_dc * ca[5];
            ca_plus_delta[2] =  ca[2] + ca[3] * delta[1] - ax1_inv * a_x_dc * ca[4];
            ca_plus_delta[4] = (ca[4] + ca[3] * delta[2] + ax1_inv * a_x_da * ca[5]) * norm_inv;
            ca_plus_delta[5] = (ca[5] + ca[3] * delta[3] - ax1_inv * a_x_da * ca[4]) * norm_inv;
        }
        else {
            const double ax1_inv = 1. / (1. - ca[3]);
            ca_plus_delta[1] =  ca[1] + ca[3] * delta[0] - ax1_inv * a_x_dc * ca[5];
            ca_plus_delta[2] =  ca[2] + ca[3] * delta[1] + ax1_inv * a_x_dc * ca[4];
            ca_plus_delta[4] = (ca[4] + ca[3] * delta[2] - ax1_inv * a_x_da * ca[5]) * norm_inv;
            ca_plus_delta[5] = (ca[5] + ca[3] * delta[3] + ax1_inv * a_x_da * ca[4]) * norm_inv;
        }

        return true;

    }

    // Jacobian of Graff13 plus operation for Ceres
    //
    virtual bool ComputeJacobian(double const* ca,
                                 double* J) const {
        
        // Jacobian has size 6x4, row-major
        // Jacobian is block diagonal: upper left and lower right 3x2 blocks are 2nd and 3rd cols of R^T, rest is zero
        J[2]  = J[3]  = J[6]  = J[7]  = J[10] = J[11] = 0.;
        J[12] = J[13] = J[16] = J[17] = J[20] = J[21] = 0.;       
        J[0] = J[14] = -ca[4];
        J[1] = J[15] = -ca[5];
        
        if (ca[3] > 0.) { // ca[3] = a[x]
            const double ax1_inv = 1. / (1. + ca[3]);
            // Jacobian has size 3x2
            J[4] = J[18] = ca[3] + ax1_inv * ca[5] * ca[5];
            J[5] = J[19] = -ax1_inv * ca[4] * ca[5];
            J[8] = J[22] = J[5];
            J[9] = J[23] = ca[3] + ax1_inv * ca[4] * ca[4];
        }
        else {
            const double ax1_inv = 1. / (1. - ca[3]);
            // Jacobian has size 3x2
            J[4] = J[18] = ca[3] - ax1_inv * ca[5] * ca[5];
            J[5] = J[19] = ax1_inv * ca[4] * ca[5];
            J[8] = J[22] = J[5];
            J[9] = J[23] = ca[3] - ax1_inv * ca[4] * ca[4];
        }
        
        return true;

    }

    virtual int GlobalSize() const { return 6; }

    virtual int LocalSize() const { return 4; }
};

#endif // LOCAL_PARAM_GRAFF13_H_
