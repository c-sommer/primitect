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

#ifndef LOCAL_PARAM_S2_H_
#define LOCAL_PARAM_S2_H_

#include <ceres/local_parameterization.h>

class LocalParamS2 : public ceres::LocalParameterization {

public:

    LocalParamS2() {}
    
    virtual ~LocalParamS2() {}

    // S2 plus operation for Ceres
    //
    virtual bool Plus(double const* n, double const* delta,
                      double* n_plus_delta) const {
        
        const double norm_inv = 1. / std::sqrt(1. + delta[0]*delta[0] + delta[1]*delta[1]);
        
        if (n[2] > 0.) {
            const double nxdelta_nz1 = (n[1] * delta[0] - n[0] * delta[1]) / (1. + n[2]);
            n_plus_delta[0] = ( nxdelta_nz1 * n[1] + n[2] * delta[0] + n[0]) * norm_inv;
            n_plus_delta[1] = (-nxdelta_nz1 * n[0] + n[2] * delta[1] + n[1]) * norm_inv;
            n_plus_delta[2] = (n[2] - n[0] * delta[0] - n[1] * delta[1]) * norm_inv;
        }
        else {
            const double nxdelta_nz1 = (n[1] * delta[0] - n[0] * delta[1]) / (1. - n[2]);
            n_plus_delta[0] = (-nxdelta_nz1 * n[1] + n[2] * delta[0] + n[0]) * norm_inv;
            n_plus_delta[1] = ( nxdelta_nz1 * n[0] + n[2] * delta[1] + n[1]) * norm_inv;
            n_plus_delta[2] = (n[2] - n[0] * delta[0] - n[1] * delta[1]) * norm_inv;
        }

        return true;

    }

    // Jacobian of S2 plus operation for Ceres
    //
    virtual bool ComputeJacobian(double const* n,
                                 double* J) const {
        
        if (n[2] > 0.) {
            const double nz1_inv = 1. / (1. + n[2]);
            // Jacobian has size 3x2
            J[0] =  n[1] * n[1] * nz1_inv + n[2];
            J[1] = -n[0] * n[1] * nz1_inv;
            J[2] =  J[1];
            J[3] =  n[0] * n[0] * nz1_inv + n[2];
            J[4] = -n[0];
            J[5] = -n[1];
        }
        else {
            const double nz1_inv = 1. / (1. - n[2]);
            // Jacobian has size 3x2
            J[0] = -n[1] * n[1] * nz1_inv + n[2];
            J[1] =  n[0] * n[1] * nz1_inv;
            J[2] =  J[1];
            J[3] = -n[0] * n[0] * nz1_inv + n[2];
            J[4] = -n[0];
            J[5] = -n[1];
        }
        
        return true;

    }

    virtual int GlobalSize() const { return 3; }

    virtual int LocalSize() const { return 2; }
};

#endif // LOCAL_PARAM_S2_H_
