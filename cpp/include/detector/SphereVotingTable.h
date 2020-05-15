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

#ifndef SPHERE_VOTING_TABLE_H_
#define SPHERE_VOTING_TABLE_H_

#include <iostream>
#include <Eigen/Dense>

template <typename T, typename VoteT>
class SphereVotingTable {

    const T Rmax_;
    const T binsize_;
    const T binsize_2_;
    const T binsize_inv_;
    
    VoteT vote_threshold_ = 10;
    
    Eigen::Matrix<VoteT, Eigen::Dynamic, 1> acc_; // accumulator array
    
    bool not_valid(T radius) {
        return (radius - binsize_2_) < 0 || radius > Rmax_ || !std::isfinite(radius);
    }
    
public:

    SphereVotingTable(T Rmax, T binsize) :
        Rmax_(Rmax),
        binsize_(binsize),
        binsize_2_(.5 * binsize_),
        binsize_inv_(1. / binsize_),
        acc_(Eigen::Matrix<VoteT, Eigen::Dynamic, 1>::Zero(static_cast<size_t>(Rmax_ * binsize_inv_ + 1.5), 1))
    {
    }
    
    void vote(T radius, bool InterpolatedWeights) { // distinguish between interpolated and integer votes
        if (not_valid(radius))
            return;
        T idx_f = radius * binsize_inv_;
        if (InterpolatedWeights) {
            size_t idx = static_cast<size_t>(idx_f + .5);
            acc_[idx-1] += (.5 + idx - idx_f);
            acc_[idx]   += (.5 + idx_f - idx);
        }
        else {
            ++acc_[static_cast<size_t>(idx_f)];
        }
    }
    
    void vote(T radius, T weight, bool InterpolatedWeights) {
        if (not_valid(radius))
            return;
        T idx_f = radius * binsize_inv_;
        if (InterpolatedWeights) {
            size_t idx = static_cast<size_t>(idx_f + .5);
            acc_[idx-1] += weight * (.5 + idx - idx_f);
            acc_[idx]   += weight * (.5 + idx_f - idx);
        }
        else {
            acc_[static_cast<size_t>(idx_f)] += weight;
        }
    }
    
    void reset() {
        acc_ *= 0;
    }
    
    void set_threshold(VoteT threshold) {
        vote_threshold_ = threshold;
    }
    
    void print_votes() const {
        for (size_t idx = 0; idx <= Rmax_ * binsize_inv_ + .5; ++idx) {
            std::cout << (idx + .5) * binsize_ << "\t" << acc_[idx] << std::endl;
        }
    }
    
    std::pair<VoteT, T> get_radius(bool MeanShift = false) const {
        int max_idx, tmp;
        VoteT max = acc_.maxCoeff(&max_idx, &tmp);
        if (max < .5 * vote_threshold_) {
            return std::make_pair(0, 0);
        }
        if (MeanShift) {
            T ms_enum = 0, ms_denom = max;
            if (max_idx > 0) {
                ms_enum  -= acc_[max_idx-1];
                ms_denom += acc_[max_idx-1];
            }
            if (max_idx+1 < Rmax_ * binsize_inv_) {
                ms_enum  += acc_[max_idx+1];
                ms_denom += acc_[max_idx+1];
            }
            T ms = ms_enum / ms_denom;
            VoteT votes = ms_denom;//(ms<0) ? (max + acc_[max_idx-1]) : (max + acc_[max_idx+1]);
            if (votes < vote_threshold_) {
                return std::make_pair(0, 0);
            }
            return std::make_pair(votes, (max_idx + ms + .5) * binsize_);
        }
        else {
            return std::make_pair(max, (max_idx + .5) * binsize_);
        }
    }

};

#endif // SPHERE_VOTING_TABLE_H_
