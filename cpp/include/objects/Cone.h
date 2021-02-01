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

#ifndef CONE_H_
#define CONE_H_

#include <iostream>
#include <Eigen/Dense>
#include "MyObject.h"

template <typename T> class Cone;
template <typename T> std::ostream& operator<<(std::ostream& os, const Cone<T>& C);

/*
 * class Cone
 */
template <typename T>
class Cone : public MyObject<T> {

    using Vec3 = Eigen::Matrix<T, 3, 1>;

    union {
        struct {
            Vec3 c_; // cone apex (3D coordinates)
            Vec3 a_; // cone axis (3D vector with unit norm)
            T sin_theta_;
            T cos_theta_;
            T theta_; // cone opening angle (scalar, in radians)
        };
        T data_[9];
    };
    
    T angle_point_axis(Vec3 a) {
        T sin_theta = a.dot(a_ - cos_theta_ * (this->rep_ - c_).normalized());
        return std::asin(sin_theta);
    }
    
    T shift_axis(Vec3 a, Vec3 c) {
        Vec3 d = this->rep_ - c;
        T h = d.dot(a);
        T r = std::sqrt(d.squaredNorm() - h * h);
        return r * cos_theta_ / sin_theta_ - h;
    }
    
public:

    Cone(Vec3 c, Vec3 a, T sin_theta, Vec3 p = Vec3(0., 0., 0.)) : // TODO: replace default rep_
        MyObject<T>(p),
        c_(c),
        a_(a.normalized()),
        sin_theta_(sin_theta),
        cos_theta_(std::sqrt(1. - sin_theta * sin_theta)),
        theta_(std::asin(sin_theta_))
    {}

    T sdf(Vec3 p) const {
        const Vec3 d = (p - c_);
        const T h = d.dot(a_);
        const T r = std::sqrt(d.squaredNorm() - h * h);
        if (h * cos_theta_ + r * sin_theta_ < 0) // point above cone apex
            return -d.norm();
        return h * sin_theta_ - r * cos_theta_;
    }
    
    Vec3 normal(Vec3 p) const {
        const Vec3 d = (p - c_);
        const T h = d.dot(a_);
        const T r = std::sqrt(d.squaredNorm() - h * h);
        if (h * cos_theta_ + r * sin_theta_ < 0) // point above cone apex
            return d.normalized();
        return -sin_theta_ * a_ + cos_theta_ * (d - h * a_).normalized();
    }
    
    Vec3 project(Vec3 p) const { // TODO : possibly improve
        const Vec3 d = (p - c_);
        const T h = d.dot(a_);
        const T r = std::sqrt(d.squaredNorm() - h * h);
        if (h * cos_theta_ + r * sin_theta_ < 0) // point above cone apex
            return c_;
        return p + sdf(p) * normal(p);
    }
    
    virtual Vec3 normal_rep() const {
        const Vec3 d = (this->rep_ - c_);
        return (-a_ + cos_theta_ * d.normalized()) / sin_theta_;
    }
    
    T integrate(T w_this, Cone* other, T w_other) {
        T w_new = w_this + w_other;
        T w_new_inv = 1. / w_new;
        c_ = (w_this * c_ + w_other * other->c_) * w_new_inv;
        a_ = (w_this * a_ + w_other * other->a_).normalized();
        theta_ = (w_this * theta_ + w_other * other->theta_) * w_new_inv;
        sin_theta_ = std::sin(theta_);
        cos_theta_ = std::cos(theta_);
        this->rep_ = project(this->rep_);
        return w_new;
    }
    
    T dist(Cone<T>* other) {
        return (c_ - other->c_).norm();
    }
    
    T angle(Cone<T>* other) {
        return a_.dot(other->a_);
    }
    
    T angle_dist(Cone<T>* other) {
        return theta_ - other->theta_;
    }

    T* data() {
        return data_;
    }
    
    template <typename U>
    Cone<U> cast() {
        return Cone<U>(c_.cast<U>(), a_.cast<U>(), U(sin_theta_), this->rep_.cast<U>());
    }
    
    friend std::ostream& operator<< <>(std::ostream& os, const Cone<T>& C);

};

/**
 * ostream definition
 */
template <typename T>
inline std::ostream& operator<<(std::ostream& os, const Cone<T>& C) {
    // display: cone apex (3D), cone axis (3D), cone opening angle (1D)
    os << C.c_.transpose() << "\t" << C.a_.transpose() << "\t" << C.theta_;
    return os;
    
}

#endif // CONE_H_
