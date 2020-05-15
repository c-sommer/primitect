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

#ifndef PLANE_CYLINDER_DETECTOR_H_
#define PLANE_CYLINDER_DETECTOR_H_

#include <iostream>
#include <fstream>
#include <map>
#include <Eigen/Dense>
#include "PpfDetector.h"
#include "CylinderVotingTable.h"
#include "objects/Plane.h"
#include "objects/Cylinder.h"

template <class Pcd, typename T, typename VoteT>
class PlaneCylinderDetector : public PpfDetector<Pcd, T, VoteT> {

    using typename PpfDetector<Pcd, T, VoteT>::Vec3;
    using typename PpfDetector<Pcd, T, VoteT>::Vec4;

    static T radius(Vec4 cppf) {
        T denom = 2. * (cppf[3] - 1);
        return (cppf[1] - cppf[2]) / denom;
    }
    
    T angle(Vec3 ni) {
        Vec3 m = this->nr_.cross(ni);
        Eigen::Matrix<T, 2, 1> m_rot = this->rot_ * m; // in this setting, no normalization necessary
        return std::atan(m_rot[1] / m_rot[0]);
    }
    
    VoteT plane_counter = 0;
    
    CylinderVotingTable<T, VoteT> cvt;
    
    std::multimap<T, Plane<T>*, std::greater<T>> plane_candidates_;
    std::multimap<T, Cylinder<T>*, std::greater<T>> cyl_candidates_;

    virtual void pair_vote(Vec4 cppf, Vec3 tmp, Vec3 ni) { // conditions: normals not parallel, feature symmetric, vectors form triangle
        DetectorSettings& s = this->settings_;
        if (cppf[3] >= this->cos_thresh_) { // if normals are approximately parallel, no Cylinder voting
            if (cppf[0] < this->min_dist_sq_) // distance between points to small to distinguish plane from cylinder/cone
                return;
            if (std::abs(cppf[1]) > this->dist_bin_ || std::abs(cppf[2]) > this->dist_bin_) // points too far from other plane
                return;
            // this already includes |angle(n, d) - pi/2| < angle_bin_
            // if sin(angle_bin_) * 10 > 1
            // which is true if angle_bin_ > 5.74° (so we require angle_bin_ >= 6° to be on the safe side)
            if (!s.ClosenessWeights) {
                if (cppf[3] >= this->cos_thresh_half_)
                    ++plane_counter;
            }
            else {
                T bin_inv = this->angle_bin_inv_;
                T d_inv = 1. / std::sqrt(cppf[0]);
                T weight = (1. - std::acos(cppf[3]) * bin_inv);
                weight *= (1. - std::abs(std::asin(cppf[1] * d_inv)) * bin_inv);
                weight *= (1. - std::abs(std::asin(cppf[2] * d_inv)) * bin_inv);
                plane_counter += weight;
            }
            return;
        }
        if (cppf[1] > 0 || cppf[2] < 0) // "quadrant/octant conditions" - possibly extend by sin_thresh or so?
            return;
        if (s.ApproximatedConditions) { // ApproximatedConditions => !PpfConditions, !ClosenessWeights
            // approximated conditions (non-uniform, non-interpretable)
            if (std::abs(cppf[1] + cppf[2]) < this->sin_thresh_) {
                cvt.vote(radius(cppf), angle(ni), s.InterpolatedWeights);
            }
        }
        else if (!s.PpfConditions) { // !ApproximatedConditions, !PpfConditions
            T S1 = std::sqrt(cppf[0] - cppf[1] * cppf[1]);
            T S2 = std::sqrt(cppf[0] - cppf[2] * cppf[2]);
            T c1 = S1 * S2 - cppf[1] * cppf[2];
            if (c1 > this->cos_thresh2_ * cppf[0]) {
                if (!s.ClosenessWeights) {
                    if (c1 > this->cos_thresh_ * cppf[0])
                        cvt.vote(radius(cppf), angle(ni), s.InterpolatedWeights);
                }
                else if (c1 < cppf[0]) { // TODO: possibly move condition further up
                    T bin_inv = this->angle_bin_inv_;
                    T weight = (1. - .5 * std::acos(c1 / cppf[0]) * bin_inv);
                    cvt.vote(radius(cppf), angle(ni), weight, s.InterpolatedWeights);
                }
            }
        }
        else {
            Vec4 ppf = PpfDetector<Pcd, T, VoteT>::compute_ppf(cppf);
            T w1 = this->angle_bin_ - .5 * std::abs(M_PI - ppf[1] - ppf[2]);
            if (w1 > 0) {
                if (!s.ClosenessWeights) { // !ApproximatedConditions, PpfConditions, !ClosenessWeights
                    cvt.vote(radius(cppf), angle(ni), s.InterpolatedWeights);
                }
                else { // !ApproximatedConditions, PpfConditions, ClosenessWeights
                    T bin_inv = this->angle_bin_inv_;
                    T weight = w1 * bin_inv;
                    cvt.vote(radius(cppf), angle(ni), weight, s.InterpolatedWeights);
                }
            }
        }
    }
    
    virtual void add_candidate(Vec3 pr, Vec3 nr) {
        using Vec2 = Eigen::Matrix<T, 2, 1>;
        std::pair<VoteT, Vec2> el = cvt.get_radius_angle(this->settings_.MeanShift);
        VoteT cyl_votes = el.first;
        if (cyl_votes > plane_counter) {
            T radius = el.second[0];
            T angle  = el.second[1];
            Vec3 c = pr - radius * nr;
            Vec3 a = this->rot_.transpose() * Vec2(std::cos(angle), std::sin(angle));
            cyl_candidates_.emplace(cyl_votes, new Cylinder<T>(c, a, radius, pr));
        }
        else if (plane_counter > this->min_votes_) {
            plane_candidates_.emplace(plane_counter, new Plane<T>(nr, -nr.dot(pr), pr));
        }
    }

public:

    PlaneCylinderDetector(DetectorSettings settings = DetectorSettings(false, false, false, false, false, false), T Rmax = 1.) :
        PpfDetector<Pcd, T, VoteT>(settings),
        cvt(CylinderVotingTable<T, VoteT>(Rmax, this->dist_bin_, this->angle_bin_))
    {
        cvt.set_threshold(this->min_votes_);
    }

    PlaneCylinderDetector(T diameter, bool is_scene, DetectorSettings settings = DetectorSettings(false, false, false, false, false, false)) :
        PpfDetector<Pcd, T, VoteT>(settings, diameter, is_scene),
        cvt(CylinderVotingTable<T, VoteT>(this->max_param_, this->dist_bin_, this->angle_bin_))
    {
        cvt.set_threshold(this->min_votes_);
    }

    virtual ~PlaneCylinderDetector() {
        for (auto& el : plane_candidates_) {
            delete el.second;
        }
        for (auto& el : cyl_candidates_) {
            delete el.second;
        }
    }
    
    virtual void cluster_candidates() {
        this->template cluster_candidates_base<Plane<T>>(plane_candidates_);
        this->template cluster_candidates_base<Cylinder<T>>(cyl_candidates_);
    }
    
    virtual void reset() {
        plane_counter = 0;
        cvt.reset();
    }
    
    virtual std::vector<MyObject<T>*> candidate_vector(const T min_vote = 0, const T max_points = std::numeric_limits<T>::infinity()) const {
    
        std::multimap<T, MyObject<T>*, std::greater<T>> candidate_map;
        
        for (auto el : plane_candidates_) {
            candidate_map.emplace(std::sqrt(el.first), static_cast<MyObject<T>*>(el.second));
        }
        for (auto el : cyl_candidates_) {
            candidate_map.emplace(std::sqrt(el.first), static_cast<MyObject<T>*>(el.second));
        }
        
        std::vector<MyObject<T>*> candidates;
        T sum_votes = 0;
        
        for (auto el : candidate_map) {
            if (el.first > min_vote && sum_votes < max_points) {
                candidates.push_back(el.second);
                sum_votes += el.first;
            }
            else
                break;
        }
        
        return candidates;
    }
    
    virtual void print_info() const {
        std::cout << plane_candidates_.size() << " plane candidates," << std::endl
                  << cyl_candidates_.size() << " cylinder candidates." << std::endl;
    }

    virtual void write_candidates(std::ostream& os) const {
        os << "# planes:" << std::endl;
        for (const auto el : plane_candidates_) {
            os << el.first << "\t" << *(el.second) << std::endl;
        }
        os << "# cylinders:" << std::endl;
        for (const auto el : cyl_candidates_) {
            os << el.first << "\t" << *(el.second) << std::endl;
        }
    }
    
    virtual void write_results(const std::string folder, const std::string filename, const T min_vote = 0, const T max_points = std::numeric_limits<T>::infinity()) const {
    
        std::string prefix = filename.substr(0, filename.size() - 4);
        
        std::ofstream plane_file(folder + prefix + "_planes.txt");
        std::ofstream sphere_file(folder + prefix + "_spheres.txt");
        sphere_file.close();
        std::ofstream cyl_file(folder + prefix + "_cylinders.txt");
        std::ofstream cone_file(folder + prefix + "_cones.txt");
        cone_file.close();
        
        std::vector<MyObject<T>*> candidates = candidate_vector(min_vote, max_points); // add second argument!!
        size_t instance = 1;
        for (const auto& obj : candidates) {        
            if (Plane<T>* p = dynamic_cast<Plane<T>*>(obj)) {
                plane_file << instance << "\t" << *p << std::endl;
            }
            else if (Cylinder<T>* c = dynamic_cast<Cylinder<T>*>(obj)) {
                cyl_file << instance << "\t" << *c << std::endl;
            }
            ++instance;
        }
        
        plane_file.close();
        cyl_file.close();
        
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

#endif // PLANE_CYLINDER_DETECTOR_H_
