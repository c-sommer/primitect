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

#ifndef PLANE_DETECTOR_H_
#define PLANE_DETECTOR_H_

#include <iostream>
#include <fstream>
#include <map>
#include <Eigen/Dense>
#include "PpfDetector.h"
#include "objects/Plane.h"

template <class Pcd, typename T, typename VoteT>
class PlaneDetector : public PpfDetector<Pcd, T, VoteT> {

    using typename PpfDetector<Pcd, T, VoteT>::Vec3;
    using typename PpfDetector<Pcd, T, VoteT>::Vec4;

    VoteT counter = 0;
    
    std::multimap<T, Plane<T>*, std::greater<T>> candidates_;

    virtual void pair_vote(Vec4 cppf, Vec3 tmp1, Vec3 tmp2) { // conditions: normals parallel, orthogonal to distance vector
        if (cppf[3] < this->cos_thresh_) // if normals are not parallel, quit
            return;
        if (cppf[0] < this->min_dist_sq_) // distance between points to small to distinguish plane from cylinder/cone
            return;
        if (std::abs(cppf[1]) > this->dist_bin_ || std::abs(cppf[2]) > this->dist_bin_) // points too far from other plane
            return;
        // this already includes |angle(n, d) - pi/2| < angle_bin_
        // if sin(angle_bin_) * 10 > 1
        // which is true if angle_bin_ > 5.74° (so we require angle_bin_ >= 6° to be on the safe side)
            
        DetectorSettings& s = this->settings_;

        if (!s.ClosenessWeights) {
            if (cppf[3] >= this->cos_thresh_half_)
                ++counter;
        }
        else {
            T bin_inv = this->angle_bin_inv_;
            T d_inv = 1. / std::sqrt(cppf[0]);
            T weight = (1. - std::acos(cppf[3]) * bin_inv);
            weight *= (1. - std::abs(std::asin(cppf[1] * d_inv)) * bin_inv);
            weight *= (1. - std::abs(std::asin(cppf[2] * d_inv)) * bin_inv);
            counter += weight;
        }
    }
    
    virtual void add_candidate(Vec3 pr, Vec3 nr) {
        if (counter > this->min_votes_) {
            candidates_.emplace(counter, new Plane<T>(nr, -nr.dot(pr), pr));
        }
    }

public:

    PlaneDetector(DetectorSettings settings = DetectorSettings(false, false, false, false, false, false)) :
        PpfDetector<Pcd, T, VoteT>(settings)
    {}
    
    PlaneDetector(T diameter, bool is_scene, DetectorSettings settings = DetectorSettings(false, false, false, false, false, false)) :
        PpfDetector<Pcd, T, VoteT>(settings, diameter, is_scene)
    {}

    virtual ~PlaneDetector() {
        for (auto& el : candidates_) {
            delete el.second;
        }
    }
    
    virtual void cluster_candidates() {
        this->template cluster_candidates_base<Plane<T>>(candidates_);
    }
    
    virtual void reset() {
        counter = 0;
    }
    
    const std::multimap<T, Plane<T>*, std::greater<T>>& candidates() const {
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
        std::cout << "Votes for Plane detector:" << counter << std::endl;
    }
    
    virtual void print_parameters(Vec3 pr, Vec3 nr) const {
        if (counter > 0)
            std::cout << "n = (" << nr.transpose() << "),\t d = " << -nr.dot(pr) << std::endl;
    }
    
    virtual void write_candidates(std::ostream& os) const {
        for (const auto el : candidates_) {
            os << el.first << "\t" << *(el.second) << std::endl;
        }
    }
    
    virtual void write_results(const std::string folder, const std::string filename, const T min_vote = 0, const T max_points = std::numeric_limits<T>::infinity()) const {
    
        std::string prefix = filename.substr(0, filename.size() - 4);
        
        std::ofstream plane_file(folder + prefix + "_planes.txt");
        std::ofstream sphere_file(folder + prefix + "_spheres.txt");
        sphere_file.close();
        std::ofstream cyl_file(folder + prefix + "_cylinders.txt");
        cyl_file.close();
        std::ofstream cone_file(folder + prefix + "_cones.txt");
        cone_file.close();
        
        std::vector<MyObject<T>*> candidates = candidate_vector(min_vote, max_points); // add second argument!!
        size_t instance = 1;
        for (const auto& obj : candidates) {        
            Plane<T>* p = dynamic_cast<Plane<T>*>(obj);
            plane_file << instance << "\t" << *p << std::endl;
            ++instance;
        }
        
        plane_file.close();
        
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

#endif // PLANE_DETECTOR_H_
