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

#ifndef PLANE_H_
#define PLANE_H_

#include <iostream>
#include <Eigen/Dense>
#include "MyObject.h"

template <typename T> class Plane;
template <typename T> std::ostream& operator<<(std::ostream& os, const Plane<T>& P);

/*
 * class Plane
 */
template <typename T>
class Plane : public MyObject<T> {

    using Vec3 = Eigen::Matrix<T, 3, 1>;

    union {
        struct {
            Vec3 n_;
            T d_;
        };
        T data_[4];
    };
    
    static int sgn(T x) {
        return (0. < x) - (x < 0.);
    }
    
public:

    Plane(Vec3 n, T d, Vec3 p = Vec3(0., 0., 0.)) : // TODO: replace default rep_
        MyObject<T>(p),
        n_(n.normalized()),
        d_(d)
    {}

    T sdf(Vec3 p) const {
        return n_.dot(p) + d_;
    }
    
    Vec3 normal(Vec3 p) const {
        return -n_;
    }
    
    Vec3 project(Vec3 p) const {
        return p - sdf(p) * n_;
    }
    
    virtual Vec3 normal_rep() const {
        return -n_;
    }
    
    T integrate(T w_this, Plane* other, T w_other) {
        // TODO: make a better plane averaging by points
        T w_new = w_this + w_other;
        T w_new_inv = 1. / w_new;
        int sign = sgn(n_.dot(other->n_));
        n_ = (w_this * n_ + sign * w_other * other->n_).normalized();
        d_ = (w_this * d_ + sign * w_other * other->d_) * w_new_inv;
        this->rep_ = project(this->rep_);
        return w_new;
    }
    
    T dist(Plane<T>* other) {
        return std::abs(d_ - other->d_);
    }
    
    T angle(Plane<T>* other) {
        return n_.dot(other->n_);
    }

    T* data() {
        return data_;
    }
    
    template <typename U>
    Plane<U> cast() {
        return Plane<U>(n_.cast<U>(), U(d_), this->rep_.cast<U>());
    }
    
    friend std::ostream& operator<< <>(std::ostream& os, const Plane<T>& P);

};

/**
 * ostream definition
 */
template <typename T>
inline std::ostream& operator<<(std::ostream& os, const Plane<T>& P) {

    os << P.n_.transpose() << "\t" << P.d_;
    return os;
    
}

#endif // PLANE_H_
