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

#ifndef CONE_DETECTOR_H_
#define CONE_DETECTOR_H_

#include <iostream>
#include <fstream>
#include <map>
#include <Eigen/Dense>
#include "PpfDetector.h"
#include "ConeVotingTable.h"
#include "objects/Cone.h"

template <class Pcd, typename T, typename VoteT>
class ConeDetector : public PpfDetector<Pcd, T, VoteT> {

    using typename PpfDetector<Pcd, T, VoteT>::Vec3;
    using typename PpfDetector<Pcd, T, VoteT>::Vec4;
    
    static T axis_dist(Vec4 cppf) {
        return cppf[2] / (1. - cppf[3]);
    }
    
    Vec3 axis(Vec4 cppf, Vec3 pi, Vec3 ni) {
        Vec3 a = this->pr_ - pi + (cppf[2] * this->nr_ + cppf[1] * ni) / (cppf[3] - 1.);
        return a.normalized();
    }
    
    ConeVotingTable<T, VoteT> cvt;
    
    std::multimap<T, Cone<T>*, std::greater<T>> candidates_;

    virtual void pair_vote(Vec4 cppf, Vec3 pi, Vec3 ni) { // conditions: normals not parallel, feature symmetric, vectors form triangle
        if (cppf[1] > 0 || cppf[2] < 0) // "quadrant/octant conditions" - possibly extend by sin_thresh or so?
            return;
        if (cppf[3] >= this->cos_thresh_) // if normals are approximately parallel, no cone voting
            return;
        DetectorSettings& s = this->settings_;
        T dist = axis_dist(cppf);
//        Vec3 a = axis(cppf, pi, ni);
//        T sin_t = std::abs(a.dot(this->nr_));
//        dist = (1. - sin_t * sin_t) * dist / sin_t;
        if (cvt.not_valid(dist))
            return;
        // do not impose any anti-degeneracy conditions
        cvt.vote(dist, axis(cppf, pi, ni), s.InterpolatedWeights);
    }
    
    virtual void add_candidate(Vec3 pr, Vec3 nr) {
        using Vec2 = Eigen::Matrix<T, 2, 1>;
        std::pair<VoteT, Vec4> el = cvt.get_dist_axis(this->settings_.MeanShift);
        VoteT votes = el.first;
        if (votes > 0) {
            T dist = el.second[0];
            Vec3 a(el.second[1], el.second[2], el.second[3]);
            Vec3 c = pr + dist * (a / a.dot(nr) - nr);
            a = (a.dot(pr - c) > 0) ? a : -a;
            T sin_theta = -a.dot(nr);
//            T sin_theta = a.dot(nr);
//            a = sin_theta > 0 ? -a : a;
//            sin_theta = std::abs(sin_theta);
//            Vec3 c = pr - dist * (a + sin_theta * nr) / (1 - sin_theta * sin_theta);
            // discard if opening angle too small (i.e. close to cylinder), or too large (i.e. close to plane)
            // but should not actually be necessary since these options are already excluded --> TODO: investigate/debug!
            if (sin_theta < this->sin_thresh_ || sin_theta > this->cos_thresh2_)
                return;
            candidates_.emplace(votes, new Cone<T>(c, a, sin_theta, pr));
        }
    }

public:

    ConeDetector(DetectorSettings settings = DetectorSettings(false, false, false, false, false, false), T Dmax = 1.) :
        PpfDetector<Pcd, T, VoteT>(settings),
        cvt(ConeVotingTable<T, VoteT>(Dmax, this->dist_bin_, this->angle_bin_))
    {
        cvt.set_threshold(this->min_votes_);
    }

    ConeDetector(T diameter, bool is_scene, DetectorSettings settings = DetectorSettings(false, false, false, false, false, false)) :
        PpfDetector<Pcd, T, VoteT>(settings, diameter, is_scene),
        cvt(ConeVotingTable<T, VoteT>(this->max_param_, this->dist_bin_, this->angle_bin_))
    {
        cvt.set_threshold(this->min_votes_);
    }

    virtual ~ConeDetector() {
        for (auto& el : candidates_) {
            delete el.second;
        }
    }
    
    virtual void cluster_candidates() {
        this->template cluster_candidates_base<Cone<T>>(candidates_);
    }
    
    virtual void reset() {
        cvt.reset();
    }
    
    const std::multimap<T, Cone<T>*, std::greater<T>>& candidates() const {
        return candidates_;
    }
    
    virtual std::vector<MyObject<T>*> candidate_vector(const T min_vote = 0, const T max_points = std::numeric_limits<T>::infinity()) const {
        std::vector<MyObject<T>*> candidates;
        T sum_votes = 0;
        T min_vote_sq = min_vote * min_vote;
        for (auto el : candidates_) {
            if (el.first > min_vote_sq && sum_votes < max_points) {
                candidates.push_back(static_cast<MyObject<T>*>(el.second));
                sum_votes += std::sqrt(el.first);
            }
            else
                break;
        }
        return candidates;
    }
    
    virtual void print_info() const {
        std::cout << candidates_.size() << " candidates." << std::endl;
    }
    
    virtual void print_votes() const {
        std::cout << "Votes for Cone detector:" << std::endl;
        cvt.print_votes();
    }
    
    virtual void write_candidates(std::ostream& os) const {
        for (const auto el : candidates_) {
            os << el.first << "\t" << *(el.second) << std::endl;
        }
    }
    
    virtual void write_results(const std::string folder, const std::string filename, const T min_vote = 0, const T max_points = std::numeric_limits<T>::infinity()) const {
    
        std::string prefix = filename.substr(0, filename.size() - 4);
        
        std::ofstream plane_file(folder + prefix + "_planes.txt");
        plane_file.close();
        std::ofstream sphere_file(folder + prefix + "_spheres.txt");
        sphere_file.close();
        std::ofstream cyl_file(folder + prefix + "_cylinders.txt");
        cyl_file.close();
        std::ofstream cone_file(folder + prefix + "_cones.txt");
        
        std::vector<MyObject<T>*> candidates = candidate_vector(min_vote, max_points); // add second argument!!
        size_t instance = 1;
        for (const auto& obj : candidates) {        
            Cone<T>* c = dynamic_cast<Cone<T>*>(obj);
            cone_file << instance << "\t" << *c << std::endl;
            ++instance;
        }
        
        cone_file.close();
        
        // write settings to file
        std::ofstream outfile(folder + prefix + "_settings.txt");

        const DetectorSettings& s = this->settings_;
        outfile << "\tInterpolatedWeights:\t"     << s.InterpolatedWeights      << std::endl
                << "\tClosenessWeights:\t"        << s.ClosenessWeights         << std::endl
                << "\tMeanShift:\t"               << s.MeanShift                << std::endl
                << "\tPpfConditions:\t"           << s.PpfConditions            << std::endl
                << "\tApproximatedConditions:\t"  << s.ApproximatedConditions   << std::endl
                << "\tClusterAveraging:\t"        << s.ClusterAveraging         << std::endl;
        outfile.close();
    }

};

#endif // CONE_DETECTOR_H_
