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


#include <ceres/ceres.h>

class TruncatedLoss : public ceres::LossFunction {

    const double scale_sq_;

public:

    TruncatedLoss(double scale = 1.) :
        scale_sq_(scale * scale)
    {}

    virtual void Evaluate(double res_sq, double out[3]) const {
        if (res_sq < scale_sq_) {
            out[0] = res_sq;
            out[1] = 1.;
            out[2] = 0.;
        }
        else {
            out[0] = scale_sq_;
            out[1] = 0.;
            out[2] = 0.;
        }
        return;
    }   
};
