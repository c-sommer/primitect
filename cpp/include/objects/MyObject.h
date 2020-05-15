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

#ifndef MY_OBJECT_H_
#define MY_OBJECT_H_

#include <iostream>
#include <Eigen/Dense>

/*
 * class Sphere
 */
template <typename T>
class MyObject {

    using Vec3 = Eigen::Matrix<T, 3, 1>;

protected:

    Vec3 rep_;
    
public:

    MyObject(Vec3 rep) :
        rep_(rep)
    {}

    virtual T sdf(Vec3 p) const = 0;
    
    virtual Vec3 normal(Vec3 p) const = 0;
    
    virtual Vec3 project(Vec3 p) const = 0;
    
    virtual Vec3 normal_rep() const = 0;
    
    virtual T* data() = 0;
    
    Vec3 rep() {
        return rep_;
    }
    
    virtual bool are_similar(MyObject<T>* other, T threshold, T cos_max_angle = std::cos(30. * M_PI / 180.)) {
        bool pt_sim = std::abs(sdf(other->rep_)) < threshold && std::abs(other->sdf(this->rep_)) < threshold;
        bool n_sim = normal(other->rep_).dot(other->normal_rep()) > cos_max_angle && other->normal(this->rep_).dot(normal_rep()) > cos_max_angle;
        return pt_sim && n_sim;
    }

};

#endif // MY_OBJECT_H_
