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

#ifndef CYLINDER_H_
#define CYLINDER_H_

#include <iostream>
#include <Eigen/Dense>
#include "MyObject.h"

template <typename T> class Cylinder;
template <typename T> std::ostream& operator<<(std::ostream& os, const Cylinder<T>& C);

/*
 * class Cylinder
 */
template <typename T>
class Cylinder : public MyObject<T> {

    using Vec3 = Eigen::Matrix<T, 3, 1>;

    union {
        struct {
            Vec3 c_; // point on cylinder axis (3D coordinates, not uniquely defined)
            Vec3 a_; // cylinder axis (3D vector with unit norm)
            T r_; // cylinder radius (scalar)
        };
        T data_[7];
    };

    static int sgn(T x) {
        return (0. < x) - (x < 0.);
    }

public:

    Cylinder(Vec3 c, Vec3 a, T r, Vec3 p = Vec3(0., 0., 0.)) : // TODO: replace default rep_
        MyObject<T>(p),
        c_(c),
        a_(a.normalized()),
        r_(r)
    {
    }

    Cylinder(const Cylinder& other) : // TODO: replace default rep_
        MyObject<T>(other.rep_),
        c_(other.c_),
        a_(other.a_),
        r_(other.r_)
    {
    }

    T sdf(Vec3 p) const {
        const Vec3 d = p - c_;
        return r_ - std::sqrt(d.squaredNorm() - d.dot(a_) * d.dot(a_));
    }

    Vec3 normal(Vec3 p) const {
        const Vec3 d = p - c_;
        return (d - a_.dot(d) * a_).normalized();
    }

    Vec3 project(Vec3 p) const {
        const Vec3 d = p - c_;
        const T dTa = d.dot(a_);
        const T sqr = d.squaredNorm() - dTa * dTa;
        return p + (r_ / std::sqrt(sqr) - 1.) * (d - dTa * a_);
    }

    virtual Vec3 normal_rep() const {
        const Vec3 d = this->rep_ - c_;
        return (d - a_.dot(d) * a_) / r_;
    }

    T integrate(T w_this, Cylinder<T>* other, T w_other) {
        T w_new = w_this + w_other;
        T w_new_inv = 1. / w_new;
        int sign = sgn(a_.dot(other->a_));
        a_ = (w_this * a_ + sign * w_other * other->a_).normalized();
        c_ = (w_this * c_ + w_other * other->c_) * w_new_inv;
        r_ = (w_this * r_ + w_other * other->r_) * w_new_inv;
        this->rep_ = project(this->rep_);
        return w_new;
    }

    T dist(Cylinder<T>* other) {
        Vec3 dc = c_ - other->c_;
        return .5 * (std::sqrt(dc.squaredNorm() - dc.dot(a_) * dc.dot(a_)) + std::sqrt(dc.squaredNorm() - dc.dot(other->a_) * dc.dot(other->a_)));
    }

    T angle(Cylinder<T>* other) {
        return a_.dot(other->a_);
    }

    T r_dist(Cylinder<T>* other) {
        return r_ - other->r_;
    }

    T* data() {
        return data_;
    }

    template <typename U>
    Cylinder<U> cast() {
        return Cylinder<U>(c_.template cast<U>(), a_.template cast<U>(), U(r_), this->rep_.template cast<U>());
    }

    friend std::ostream& operator<< <>(std::ostream& os, const Cylinder<T>& C);

};

/**
 * ostream definition
 */
template <typename T>
inline std::ostream& operator<<(std::ostream& os, const Cylinder<T>& C) {
    // display: point on cylinder axis (3D), cylinder axis (3D), cylinder radius (1D)
    os << C.c_.transpose() << "\t" << C.a_.transpose() << "\t" << C.r_;
    return os;

}

#endif // CYLINDER_H_
