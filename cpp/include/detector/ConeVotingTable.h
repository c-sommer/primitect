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

#ifndef CONE_VOTING_TABLE_H_
#define CONE_VOTING_TABLE_H_

#include <iostream>
#include <Eigen/Dense>

template <typename T>
Eigen::Matrix<size_t, 1, Eigen::Dynamic> sort_indices(const Eigen::Matrix<T, Eigen::Dynamic, 1>& vec) {

  // initialize original index locations
  Eigen::Matrix<size_t, 1, Eigen::Dynamic> idx = Eigen::Matrix<size_t, Eigen::Dynamic, 1>::Zero(vec.size());
  std::iota(idx.data(), idx.data() + idx.size(), 0);

  // sort indices based on comparing values in v
  std::partial_sort(idx.data(), idx.data() + 5, idx.data() + idx.size(),
       [&vec](size_t i1, size_t i2) { return vec[i1] > vec[i2]; });

  return idx;
}

template <typename T, typename VoteT>
class ConeVotingTable {

    const T Dmax_;
    const T binsize_;
    const T binsize_2_;
    const T binsize_inv_;

    const size_t num_angles_;
    size_t num_tot_angles_;
    
    const T anglebin_;
    const T anglebin_inv_;
    
    VoteT vote_threshold_ = 10;
    
    Eigen::Matrix<size_t, Eigen::Dynamic, 1> num_phi_;
    Eigen::Matrix<T, Eigen::Dynamic, 1> bin_phi_inv_;
    Eigen::Matrix<size_t, Eigen::Dynamic, 1> sum_num_phi_;
    
    Eigen::Matrix<VoteT, Eigen::Dynamic, Eigen::Dynamic> acc_; // accumulator array
    
    // TODO: array for neighbors! (8 neighbors, of which find 3 smallest) - possibly transpose??
    Eigen::Matrix<size_t, 4, Eigen::Dynamic> neighbors_;
    // TODO: array containing actual axes - possibly transpose??
    Eigen::Matrix<T, 3, Eigen::Dynamic> axes_;
    
public:

    ConeVotingTable(T Dmax, T binsize, T anglebin) :
        Dmax_(Dmax),
        binsize_(binsize),
        binsize_2_(.5 * binsize_),
        binsize_inv_(1. / binsize_),
        num_angles_(static_cast<size_t>(M_PI / anglebin + .5)),
        anglebin_(M_PI / num_angles_),
        anglebin_inv_(1. / anglebin_),
        num_phi_(Eigen::Matrix<size_t, Eigen::Dynamic, 1>::Zero(num_angles_, 1)),
        bin_phi_inv_(Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(num_angles_, 1)),
        sum_num_phi_(Eigen::Matrix<size_t, Eigen::Dynamic, 1>::Zero(num_angles_, 1))
    {
        num_tot_angles_ = 0;
        // TODO: create index arrays
        for (size_t i=0; i<num_angles_; ++i) {
            T theta = (.5 + i) * anglebin_;
            num_phi_[i] = static_cast<size_t>(std::sin(theta) * num_angles_ + .5);
            bin_phi_inv_[i] = num_phi_[i] / M_PI;
            sum_num_phi_[i] = num_tot_angles_;
            num_tot_angles_ += num_phi_[i];
        }
        // TODO: create accumulator
        acc_ = Eigen::Matrix<VoteT, Eigen::Dynamic, Eigen::Dynamic>::Zero(static_cast<size_t>(Dmax_ * binsize_inv_ + 1.5), num_tot_angles_);
        neighbors_ = Eigen::Matrix<size_t, 4, Eigen::Dynamic>::Zero(4, num_tot_angles_);
        axes_ = Eigen::Matrix<T, 3, Eigen::Dynamic>::Zero(3, num_tot_angles_);
        // fill axes
        size_t ij = 0;
        for (size_t i=0; i<num_angles_; ++i) {
            for (size_t j=0; j<num_phi_[i]; ++j) {
                T theta = (.5 + i) * anglebin_;
                T phi = (.5 + j) * M_PI / num_phi_[i];
                axes_(0, ij) = std::cos(phi) * std::sin(theta);
                axes_(1, ij) = std::sin(phi) * std::sin(theta);
                axes_(2, ij) = std::cos(theta);
                ++ij;
            }
        }
        // fill neighbors
        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> cos_dist = axes_.transpose() * axes_;
        cos_dist = cos_dist.cwiseAbs(); 
        for (size_t i=0; i<num_tot_angles_; ++i) {
            Eigen::Matrix<size_t, Eigen::Dynamic, 1> idx = sort_indices<T>(cos_dist.col(i));
            neighbors_.col(i) = idx.segment<4>(1);
        }
    }
    
    void vote(T axis_dist, Eigen::Matrix<T, 3, 1> axis, bool InterpolatedWeights) { // distinguish between interpolated and integer votes
        T idx_f = axis_dist * binsize_inv_;
        T theta = (axis[1] > 0) ? std::acos(axis[2]) : (M_PI - std::acos(axis[2]));
        T phi = std::atan(axis[1] / axis[0]);
        if (phi < 0)
            phi += M_PI;
        T idy_f = theta * anglebin_inv_;
        if (InterpolatedWeights) {
            size_t idx = static_cast<size_t>(idx_f + .5);
            T wx_minus = .5 + idx - idx_f;
            T wx_plus  = .5 + idx_f - idx;
            if (idy_f < .5) {
                T cos_phi = std::cos(phi);
                T sin_phi = std::sqrt(1. - cos_phi * cos_phi);
                acc_(idx-1, 0) += wx_minus * .25 * (1 + sin_phi) * (1 + cos_phi);
                acc_(idx-1, 1) += wx_minus * .25 * (1 + sin_phi) * (1 - cos_phi);
                acc_(idx-1, num_tot_angles_-2) += wx_minus * .25 * (1 - sin_phi) * (1 - cos_phi);
                acc_(idx-1, num_tot_angles_-1) += wx_minus * .25 * (1 - sin_phi) * (1 + cos_phi);
                acc_(idx, 0) += wx_plus * .25 * (1 + sin_phi) * (1 + cos_phi);
                acc_(idx, 1) += wx_plus * .25 * (1 + sin_phi) * (1 - cos_phi); 
                acc_(idx, num_tot_angles_-2) += wx_plus * .25 * (1 - sin_phi) * (1 - cos_phi);
                acc_(idx, num_tot_angles_-1) += wx_plus * .25 * (1 - sin_phi) * (1 + cos_phi);
            }
            else if (idy_f >= num_angles_ - .5) {
                T cos_phi = std::cos(phi);
                T sin_phi = std::sqrt(1. - cos_phi * cos_phi);
                acc_(idx-1, 0) += wx_minus * .25 * (1 - sin_phi) * (1 - cos_phi);
                acc_(idx-1, 1) += wx_minus * .25 * (1 - sin_phi) * (1 + cos_phi);
                acc_(idx-1, num_tot_angles_-2) += wx_minus * .25 * (1 + sin_phi) * (1 + cos_phi);
                acc_(idx-1, num_tot_angles_-1) += wx_minus * .25 * (1 + sin_phi) * (1 - cos_phi);
                acc_(idx, 0) += wx_plus * .25 * (1 - sin_phi) * (1 - cos_phi);
                acc_(idx, 1) += wx_plus * .25 * (1 - sin_phi) * (1 + cos_phi); 
                acc_(idx, num_tot_angles_-2) += wx_plus * .25 * (1 + sin_phi) * (1 + cos_phi);
                acc_(idx, num_tot_angles_-1) += wx_plus * .25 * (1 + sin_phi) * (1 - cos_phi);
            }
            else {
                size_t idy = static_cast<size_t>(idy_f + .5);
                T idz_f = phi * bin_phi_inv_[idy-1];
                size_t idz = static_cast<size_t>(idz_f + .5);
                T wyz_tmp = (.5 + idy - idy_f) * (.5 + idz - idz_f);
                if (idz_f > .5) {
                    acc_(idx-1, sum_num_phi_[idy-1] + idz-1) += wx_minus * wyz_tmp;
                    acc_(idx, sum_num_phi_[idy-1] + idz-1)   += wx_plus * wyz_tmp;
                }
                else {
                    // TODO: possibly make sum_num_phi_ one dim larger and combine the two summands in the second argument
                    acc_(idx-1, sum_num_phi_[num_angles_-idy] + num_phi_[num_angles_-idy]-1) += wx_minus * wyz_tmp;
                    acc_(idx, sum_num_phi_[num_angles_-idy] + num_phi_[num_angles_-idy]-1)   += wx_plus * wyz_tmp;
                }
                wyz_tmp = (.5 + idy - idy_f) * (.5 + idz_f - idz);
                if (idz_f < num_phi_[idy] - .5) {
                    acc_(idx-1, sum_num_phi_[idy-1] + idz) += wx_minus * wyz_tmp;
                    acc_(idx, sum_num_phi_[idy-1] + idz)   += wx_plus * wyz_tmp;
                }
                else {
                    acc_(idx-1, sum_num_phi_[num_angles_-idy]) += wx_minus * wyz_tmp;
                    acc_(idx, sum_num_phi_[num_angles_-idy])   += wx_plus * wyz_tmp;
                }
                idz_f = phi * bin_phi_inv_[idy];
                idz = static_cast<size_t>(idz_f + .5);
                wyz_tmp = (.5 + idy_f - idy) * (.5 + idz - idz_f);
                if (idz_f > .5) {
                    acc_(idx-1, sum_num_phi_[idy] + idz-1) += wx_minus * wyz_tmp;
                    acc_(idx, sum_num_phi_[idy] + idz-1)   += wx_plus * wyz_tmp;
                }
                else {
                    // TODO: possibly make sum_num_phi_ one dim larger and combine the two summands in the second argument
                    acc_(idx-1, sum_num_phi_[num_angles_-idy-1] + num_phi_[num_angles_-idy]-1) += wx_minus * wyz_tmp;
                    acc_(idx, sum_num_phi_[num_angles_-idy-1] + num_phi_[num_angles_-idy]-1)   += wx_plus * wyz_tmp;
                }
                wyz_tmp = (.5 + idy_f - idy) * (.5 + idz_f - idz);
                if (idz_f < num_phi_[idy] - .5) {
                    acc_(idx-1, sum_num_phi_[idy] + idz) += wx_minus * wyz_tmp;
                    acc_(idx, sum_num_phi_[idy] + idz)   += wx_plus * wyz_tmp;
                }
                else {
                    acc_(idx-1, sum_num_phi_[num_angles_-idy-1]) += wx_minus * wyz_tmp;
                    acc_(idx, sum_num_phi_[num_angles_-idy-1])   += wx_plus * wyz_tmp;
                }
            }
        }
        else {
            size_t idy = static_cast<size_t>(idy_f) % num_angles_;
            T idz_f = phi * bin_phi_inv_[idy];
            // TODO: add modulus to make safe
            ++acc_(static_cast<size_t>(idx_f), sum_num_phi_[idy] + static_cast<size_t>(idz_f));
        }
    }
    
    // No ClosenessWeights implementation needed, since no "equality constraints" given
    
    void reset() {
        acc_ *= 0;
    }
    
    void set_threshold(VoteT threshold) {
        vote_threshold_ = threshold;
    }
        
    bool not_valid(T axis_dist) {
        return (axis_dist - binsize_2_) < 0 || axis_dist > Dmax_ || !std::isfinite(axis_dist); // TODO: second option only if large-angle cones expected
    }
    
    void print_votes() const {
        std::cout << acc_ << std::endl;
    }
    
    std::pair<VoteT, Eigen::Matrix<T, 4, 1>> get_dist_axis(bool MeanShift = false) const {
        using Vec4 = Eigen::Matrix<T, 4, 1>;
        int max_idx, max_idy;
        VoteT max = acc_.maxCoeff(&max_idx, &max_idy);
        if (max < .125 * vote_threshold_) {
            return std::make_pair(0, Vec4::Zero());
        }
        if (MeanShift) {
            T msx_enum = 0, msx_denom = max;
            if (max_idx > 0) {
                msx_enum  -= acc_(max_idx-1, max_idy);
                msx_denom += acc_(max_idx-1, max_idy);
            }
            if (max_idx+1 < Dmax_ * binsize_inv_) {
                msx_enum  += acc_(max_idx+1, max_idy);
                msx_denom += acc_(max_idx+1, max_idy);
            }
            T msx = msx_enum / msx_denom;
            Eigen::Matrix<T, 3, 1> axis = max * axes_.col(max_idy);
            T msy_denom = max;
            for (size_t k=0; k<4; ++k) {
                size_t n_idy = neighbors_(k, max_idy);
                int sign = (axis.dot(axes_.col(n_idy)) > 0) ? 1 : -1;
                T weight = acc_(max_idx, n_idy);
                msy_denom += weight;
                axis += sign * weight * axes_.col(n_idy);
            }
            axis.normalize();
            VoteT votes = msx_denom + msy_denom - max; // TODO: proper mean shift !!
            if (votes < vote_threshold_) {
                return std::make_pair(0, Vec4::Zero());
            }
            return std::make_pair(votes, Vec4((max_idx + msx + .5) * binsize_, axis[0], axis[1], axis[2]));
        }
        else {
            return std::make_pair(max, Vec4((max_idx + .5) * binsize_, axes_(0, max_idy), axes_(1, max_idy), axes_(2, max_idy)));
        }
    }

};

#endif // CONE_VOTING_TABLE_H_
