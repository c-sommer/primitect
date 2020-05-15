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

#ifndef PPF_DETECTOR_H_
#define PPF_DETECTOR_H_

#include <iostream>
#include <bitset>
#include <Eigen/Dense>
#include "objects/MyObject.h"

struct DetectorSettings {

    bool InterpolatedWeights = false; // 1, must be false for integer voting types
    bool ClosenessWeights = false; // 2, must be false for integer voting types
    bool MeanShift = false; // 4
    bool PpfConditions = false; // 8
    bool ApproximatedConditions = false; // 16
    bool ClusterAveraging = false; // 32
    
    DetectorSettings(bool InterpolatedWeights,
                     bool ClosenessWeights,
                     bool MeanShift,
                     bool PpfConditions,
                     bool ApproximatedConditions,
                     bool ClusterAveraging) :
        InterpolatedWeights(InterpolatedWeights),
        ClosenessWeights(ClosenessWeights),
        MeanShift(MeanShift),
        PpfConditions(PpfConditions),
        ApproximatedConditions(ApproximatedConditions),
        ClusterAveraging(ClusterAveraging)
    {}
    
    DetectorSettings(std::bitset<6> all) :
        InterpolatedWeights(all[0]),
        ClosenessWeights(all[1]),
        MeanShift(all[2]),
        PpfConditions(all[3]),
        ApproximatedConditions(all[4]),
        ClusterAveraging(all[5])
    {}
    
    void print() {
        std::cout << "\tInterpolatedWeights:\t"     << InterpolatedWeights      << std::endl
                  << "\tClosenessWeights:\t"        << ClosenessWeights         << std::endl
                  << "\tMeanShift:\t"               << MeanShift                << std::endl
                  << "\tPpfConditions:\t"           << PpfConditions            << std::endl
                  << "\tApproximatedConditions:\t"  << ApproximatedConditions   << std::endl
                  << "\tClusterAveraging:\t"        << ClusterAveraging         << std::endl;
    }

};

template <class Pcd, typename T, typename VoteT>
class PpfDetector {

protected:

    using Vec3 = Eigen::Matrix<T, 3, 1>;
    using Vec4 = Eigen::Matrix<T, 4, 1>;
    
    const T angle_bin_; // size of a bin in angular direction (in rad)
    const T angle_bin_inv_;
    const T cos_thresh_;
    const T cos_thresh2_;
    const T cos_thresh_half_;
    const T sin_thresh_;
    const T dist_bin_;
    const T dist_bin_inv_;
    const T min_dist_sq_;
    const T dist_thresh_sq_;
    
    const T max_param_;
    
    const VoteT min_votes_ = 8; //2048, 8; 4096, 24
    
    static Vec4 compute_cppf(Vec3 pr, Vec3 nr, Vec3 pi, Vec3 ni) {
        Vec3 d = pi - pr;
        Vec4 cppf;
        cppf[0] = d.squaredNorm();
        cppf[1] = d.dot(nr);
        cppf[2] = d.dot(ni);
        cppf[3] = nr.dot(ni);
        return cppf;
    }
    
    static Vec4 compute_ppf(Vec4 cppf) {
        Vec4 ppf;
        ppf[0] = std::sqrt(cppf[0]);
        T ppf0_inv = 1. / ppf[0];
        ppf[1] = std::acos(cppf[1] * ppf0_inv);
        ppf[2] = std::acos(cppf[2] * ppf0_inv);
        ppf[3] = std::acos(cppf[3]);
        return ppf;
    }
    
    Vec3 pr_, nr_;
    
    Eigen::Matrix<T, 2, 3> rot_;
    
    virtual void pair_vote(Vec4 cppf, Vec3 pi, Vec3 ni) = 0;
    
    virtual void add_candidate(Vec3 pr, Vec3 nr) = 0;

    DetectorSettings settings_;

public:

    PpfDetector(DetectorSettings settings, T angle_bin = 10. * M_PI / 180., T dist_bin = 0.025) :
        settings_(settings),
        angle_bin_(angle_bin),
        angle_bin_inv_(1. / angle_bin),
        cos_thresh_(std::cos(angle_bin)),
        cos_thresh2_(std::cos(2 * angle_bin)),
        cos_thresh_half_(std::cos(.5 * angle_bin)),
        sin_thresh_(std::sin(angle_bin)),
        dist_bin_(dist_bin),
        dist_bin_inv_(1. / dist_bin),
        min_dist_sq_(64 * dist_bin * dist_bin),
        dist_thresh_sq_(1.0),
        max_param_(1.0)
    {
        correct_settings();
    }
    
    PpfDetector(DetectorSettings settings, T diameter, bool is_scene, T angle_bin = 10. * M_PI / 180.) :
        settings_(settings),
        angle_bin_(angle_bin),
        angle_bin_inv_(1. / angle_bin),
        cos_thresh_(std::cos(angle_bin)),
        cos_thresh2_(std::cos(2 * angle_bin)),
        cos_thresh_half_(std::cos(.5 * angle_bin)),
        sin_thresh_(std::sin(angle_bin)),
        dist_bin_(is_scene ? .005 * diameter : .01 * diameter),
        dist_bin_inv_(1. / dist_bin_),
        min_dist_sq_(.0016 * diameter * diameter),
        dist_thresh_sq_(is_scene ? .04 * diameter * diameter : .16 * diameter * diameter),
        max_param_(is_scene ? .2 * diameter : .5 * diameter)
    {
        correct_settings();
    }

    void correct_settings() {
        if (std::is_integral<T>::value) { // if voting type is integer, weighted voting does not make sense
            settings_.InterpolatedWeights = false;
            settings_.ClosenessWeights = false;
        }
        if (settings_.PpfConditions || settings_.ClosenessWeights) { // only CPPF conditions can be approximated in a meaningful way
            settings_.ApproximatedConditions = false;
        }    
    }

    virtual ~PpfDetector() {}

    virtual void cast_votes(Pcd& pc) {
        size_t num_points = pc.num_points();
        for (size_t r = 0; r < num_points; ++r) { // reference point
            if (!set_rot(pc.normal(r))) {
                continue;
            }
            pr_ = pc.point(r);
            nr_ = pc.normal(r);
            for (size_t i = 0; i < num_points; ++i) { // pair point
                if (i == r)
                    continue;
                Vec4 cppf = compute_cppf(pr_, nr_, pc.point(i), pc.normal(i));
                if (cppf[0] < dist_thresh_sq_)
                    pair_vote(cppf, pc.point(i), pc.normal(i));
            }
            add_candidate(pr_, nr_);
            reset();
        }
    }
    
    template <class Object>
    void cluster_candidates_base(std::multimap<T, Object*, std::greater<T>>& candidates) {
        static_assert(std::is_base_of<MyObject<T>, Object>::value, "Object must be derived from MyObject<T>"); // make sure Object is of type MyObject<T>
        bool ClusterAveraging = settings_.ClusterAveraging;
        std::multimap<T, Object*, std::greater<T>> new_candidates;
        for (auto it = candidates.begin(); it != candidates.end(); ++it) {
            Object* Oi = it->second;
            T wcurr = it->first;
            for (auto jt = std::next(it); jt != candidates.end();) {
                Object* Oj = jt->second;
                if (Oi->are_similar(Oj, .1)) {
                    if (ClusterAveraging) // if false: simple non-maximum suppression
                        wcurr = Oi->integrate(wcurr, Oj, jt->first);
                    else
                        wcurr += jt->first;
                    jt = candidates.erase(jt);
                } else {
                    ++jt;
                }
            }
            new_candidates.emplace(wcurr, Oi);
        }
        candidates = new_candidates;
        return;
    }
    
    virtual void cluster_candidates() = 0;
    
    virtual void reset() = 0;
    
    void set_settings(DetectorSettings settings) {
        settings_ = settings;
    }
    
    bool set_rot(Vec3 nr) {
        const T nx1_inv = 1. / (1. + nr[0]);
        if (!std::isfinite(nx1_inv)) {
            return false;
        }
        rot_(0,0) = -nr[1];
        rot_(1,0) = -nr[2];
        rot_(0,1) = 1. - nr[1] * nr[1] * nx1_inv;
        rot_(1,1) = -nr[1] * nr[2] * nx1_inv;
        rot_(0,2) = rot_(1,1);
        rot_(1,2) = 1. - nr[2] * nr[2] * nx1_inv;
        return true;
    }
    
    virtual void print_votes() const {
        std::cout << "Votes cannot be printed (yet), sorry." << std::endl;
    }
    
    virtual void print_parameters(Vec3 pr, Vec3 nr) const {
        std::cout << "This functionality is not (yet) implemented, sorry." << std::endl;
    }
    
    virtual void print_info() const = 0;
    
    virtual void write_candidates(std::ostream& os) const = 0;
    
    virtual void print_candidates(const int threshold = 0) const {
        write_candidates(std::cout);
    }
    
    virtual void candidates_to_file(std::string filename, const int threshold = 0) const {
        std::ofstream outfile(filename);
        write_candidates(outfile);
        outfile.close();
    }
    
    virtual void write_results(const std::string folder, std::string filename, const T min_vote = 0, const T max_points = std::numeric_limits<T>::infinity()) const {}
    
    virtual std::vector<MyObject<T>*> candidate_vector(const T min_vote = 0, const T max_points = std::numeric_limits<T>::infinity()) const = 0;

};

#endif // PPF_DETECTOR_H_
