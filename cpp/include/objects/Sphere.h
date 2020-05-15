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

#ifndef SPHERE_H_
#define SPHERE_H_

#include <iostream>
#include <Eigen/Dense>
#include "MyObject.h"

template <typename T> class Sphere;
template <typename T> std::ostream& operator<<(std::ostream& os, const Sphere<T>& S);

/*
 * class Sphere
 */
template <typename T>
class Sphere : public MyObject<T> {

    using Vec3 = Eigen::Matrix<T, 3, 1>;

    union {
        struct {
            Vec3 c_;
            T r_;
        };
        T data_[4];
    };
    
public:

    Sphere(Vec3 c, T r, Vec3 p = Vec3(0., 0., 0.)) : // TODO: replace default rep_
        MyObject<T>(p),
        c_(c),
        r_(r)
    {}

    T sdf(Vec3 p) const {
        const Vec3 d = p - c_;
        return r_ - d.norm();
    }
    
    Vec3 normal(Vec3 p) const {
        return (p - c_).normalized();
    }
    
    Vec3 project(Vec3 p) const {
        const Vec3 d = p - c_;
        return p + (r_ / d.norm() - 1.) * d;
    }
    
    virtual Vec3 normal_rep() const {
        return (this->rep_ - c_) / r_;
    }
    
    T integrate(T w_this, Sphere* other, T w_other) {
        T w_new = w_this + w_other;
        T w_new_inv = 1. / w_new;
        c_ = (w_this * c_ + w_other * other->c_) * w_new_inv;
        r_ = (w_this * r_ + w_other * other->r_) * w_new_inv;
        this->rep_ = project(this->rep_);
        return w_new;
    }
    
    T dist(Sphere<T>* other) {
        return (c_ - other->c_).norm();
    }
    
    T r_dist(Sphere<T>* other) {
        return r_ - other->r_;
    }

    T* data() {
        return data_;
    }
    
    template <typename U>
    Sphere<U> cast() {
        return Sphere<U>(c_.cast<U>(), U(r_), this->rep_.cast<U>());
    }
    
    friend std::ostream& operator<< <>(std::ostream& os, const Sphere<T>& S);

};

/**
 * ostream definition
 */
template <typename T>
inline std::ostream& operator<<(std::ostream& os, const Sphere<T>& S) {

    os << S.c_.transpose() << "\t" << S.r_;  
    return os;
    
}

#endif // SPHERE_H_
