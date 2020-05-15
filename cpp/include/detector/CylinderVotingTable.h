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

#ifndef CYLINDER_VOTING_TABLE_H_
#define CYLINDER_VOTING_TABLE_H_

#include <iostream>
#include <Eigen/Dense>

template <typename T, typename VoteT>
class CylinderVotingTable {

    const T Rmax_;
    const T binsize_;
    const T binsize_2_;
    const T binsize_inv_;

    const size_t num_angles_;
    
    const T anglebin_;
    const T anglebin_inv_;
    
    VoteT vote_threshold_ = 10;
    
    Eigen::Matrix<VoteT, Eigen::Dynamic, Eigen::Dynamic> acc_; // accumulator array
    
    bool not_valid(T radius) {
        return (radius - binsize_2_) < 0 || radius > Rmax_ || !std::isfinite(radius);
    }
    
public:

    CylinderVotingTable(T Rmax, T binsize, T anglebin) :
        Rmax_(Rmax),
        binsize_(binsize),
        binsize_2_(.5 * binsize_),
        binsize_inv_(1. / binsize_),
        num_angles_(static_cast<size_t>(M_PI / anglebin + .5)),
        anglebin_(M_PI / num_angles_),
        anglebin_inv_(1. / anglebin_),
        acc_(Eigen::Matrix<VoteT, Eigen::Dynamic, Eigen::Dynamic>::Zero(static_cast<size_t>(Rmax_ * binsize_inv_ + 1.5), num_angles_))
    {
    }
    
    void vote(T radius, T angle, bool InterpolatedWeights) { // distinguish between interpolated and integer votes
        if (not_valid(radius) || !std::isfinite(angle))
            return;
        T idx_f = radius * binsize_inv_;
        T idy_f = angle * anglebin_inv_;
        if (idy_f < 0)
            idy_f += M_PI * anglebin_inv_; // TODO: replace by num_angles_
        if (InterpolatedWeights) {
            size_t idx = static_cast<size_t>(idx_f + .5);
            size_t idy = static_cast<size_t>(idy_f + .5);
            if (idy > 0) {
                acc_(idx-1, idy-1) += (.5 + idx - idx_f) * (.5 + idy - idy_f);
                acc_(idx, idy-1)   += (.5 + idx_f - idx) * (.5 + idy - idy_f);
            }
            else {
                acc_(idx-1, num_angles_-1) += (.5 + idx - idx_f) * (.5 + idy - idy_f);
                acc_(idx, num_angles_-1)   += (.5 + idx_f - idx) * (.5 + idy - idy_f);
            }
            if (idy < num_angles_) {
                acc_(idx-1, idy) += (.5 + idx - idx_f) * (.5 + idy_f - idy);
                acc_(idx, idy)   += (.5 + idx_f - idx) * (.5 + idy_f - idy);
            }
            else {
                acc_(idx-1, 0) += (.5 + idx - idx_f) * (.5 + idy_f - idy);
                acc_(idx, 0)   += (.5 + idx_f - idx) * (.5 + idy_f - idy);
            }
        }
        else {
            ++acc_(static_cast<size_t>(idx_f), static_cast<size_t>(idy_f) % num_angles_);
        }
    }
    
    void vote(T radius, T angle, T weight, bool InterpolatedWeights) {
        if (not_valid(radius) || !std::isfinite(angle))
            return;
        T idx_f = radius * binsize_inv_;
        T idy_f = angle * anglebin_inv_;
        if (idy_f < 0)
            idy_f += M_PI * anglebin_inv_; // TODO: replace by num_angles_
        if (InterpolatedWeights) {
            size_t idx = static_cast<size_t>(idx_f + .5);
            size_t idy = static_cast<size_t>(idy_f + .5);
            if (idy > 0) {
                acc_(idx-1, idy-1) += weight * (.5 + idx - idx_f) * (.5 + idy - idy_f);
                acc_(idx, idy-1)   += weight * (.5 + idx_f - idx) * (.5 + idy - idy_f);
            }
            else {
                acc_(idx-1, num_angles_-1) += weight * (.5 + idx - idx_f) * (.5 + idy - idy_f);
                acc_(idx, num_angles_-1)   += weight * (.5 + idx_f - idx) * (.5 + idy - idy_f);
            }
            if (idy < num_angles_) {
                acc_(idx-1, idy) += weight * (.5 + idx - idx_f) * (.5 + idy_f - idy);
                acc_(idx, idy)   += weight * (.5 + idx_f - idx) * (.5 + idy_f - idy);
            }
            else {
                acc_(idx-1, 0) += weight * (.5 + idx - idx_f) * (.5 + idy_f - idy);
                acc_(idx, 0)   += weight * (.5 + idx_f - idx) * (.5 + idy_f - idy);
            }
        }
        else {
            acc_(static_cast<size_t>(idx_f), static_cast<size_t>(idy_f) % num_angles_) += weight;
        }
    }
    
    void reset() {
        acc_ *= 0;
    }
    
    void set_threshold(VoteT threshold) {
        vote_threshold_ = threshold;
    }
    
    void print_votes() const {
        std::cout << acc_ << std::endl;
    }
    
    std::pair<VoteT, Eigen::Matrix<T, 2, 1>> get_radius_angle(bool MeanShift = false) const {
        using Vec2 = Eigen::Matrix<T, 2, 1>;
        int max_idx, max_idy;
        VoteT max = acc_.maxCoeff(&max_idx, &max_idy);
        if (max < .25 * vote_threshold_) {
            return std::make_pair(0, Vec2(0, 0));
        }
        if (MeanShift) {
            T msx_enum = 0, msx_denom = max, msy_enum = 0, msy_denom = max;
            if (max_idx > 0) {
                msx_enum  -= acc_(max_idx-1, max_idy);
                msx_denom += acc_(max_idx-1, max_idy);
            }
            if (max_idx+1 < Rmax_ * binsize_inv_) {
                msx_enum  += acc_(max_idx+1, max_idy);
                msx_denom += acc_(max_idx+1, max_idy);
            }
            if (max_idy > 0) {
                msy_enum  -= acc_(max_idx, max_idy-1);
                msy_denom += acc_(max_idx, max_idy-1);
            } else {
                msy_enum  -= acc_(max_idx, acc_.cols()-1);
                msy_denom -= acc_(max_idx, acc_.cols()-1);
            }
            if (max_idy+1 < acc_.cols()) {
                msy_enum  += acc_(max_idx, max_idy+1);
                msy_denom += acc_(max_idx, max_idy+1);
            } else {
                msy_enum  += acc_(max_idx, 0);
                msy_denom += acc_(max_idx, 0);
            }
            T msx = msx_enum / msx_denom;
            T msy = msy_enum / msy_denom;
            VoteT votes = msx_denom + msy_denom - max; // TODO: proper mean shift !! //(ms<0) ? (max + acc_[max_idx-1]) : (max + acc_[max_idx+1]);
            if (votes < vote_threshold_) {
                return std::make_pair(0, Vec2(0, 0));
            }
            return std::make_pair(votes, Vec2((max_idx + msx + .5) * binsize_, (max_idy + msy +.5) * anglebin_));
        }
        else {
            return std::make_pair(max, Vec2((max_idx + .5) * binsize_, (max_idy + .5) * anglebin_));
        }
    }

};

#endif // CYLINDER_VOTING_TABLE_H_
