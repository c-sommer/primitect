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

#ifndef PRIMITIVE_DETECTOR_H_
#define PRIMITIVE_DETECTOR_H_

#include <iostream>
#include <fstream>
#include <map>
#include <Eigen/Dense>
#include "PpfDetector.h"
#include "CylinderVotingTable.h"
#include "objects/Plane.h"
#include "objects/Cylinder.h"

template <class Pcd, typename T, typename VoteT>
class PrimitiveDetector : public PpfDetector<Pcd, T, VoteT> {

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
    
    static T axis_dist(Vec4 cppf) {
        return cppf[2] / (1. - cppf[3]);
    }
    
    Vec3 axis(Vec4 cppf, Vec3 pi, Vec3 ni) {
        Vec3 a = this->pr_ - pi + (cppf[2] * this->nr_ + cppf[1] * ni) / (cppf[3] - 1.);
        return a.normalized();
    }
    
    VoteT plane_counter = 0;
    SphereVotingTable<T, VoteT> svt;
    CylinderVotingTable<T, VoteT> cyvt;
    ConeVotingTable<T, VoteT> covt;
    
    std::multimap<T, Plane<T>*, std::greater<T>> plane_candidates_;
    std::multimap<T, Sphere<T>*, std::greater<T>> sphere_candidates_;
    std::multimap<T, Cylinder<T>*, std::greater<T>> cyl_candidates_;
    std::multimap<T, Cone<T>*, std::greater<T>> cone_candidates_;

    virtual void pair_vote(Vec4 cppf, Vec3 pi, Vec3 ni) {
        DetectorSettings& s = this->settings_;
        // ---------- PLANE ----------
        if (cppf[3] >= this->cos_thresh_) { // if normals are approximately parallel, check rest of plane conditions
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
                const T bin_inv = this->angle_bin_inv_;
                const T d_inv = 1. / std::sqrt(cppf[0]);
                T weight = (1. - std::acos(cppf[3]) * bin_inv);
                weight *= (1. - std::abs(std::asin(cppf[1] * d_inv)) * bin_inv);
                weight *= (1. - std::abs(std::asin(cppf[2] * d_inv)) * bin_inv);
                plane_counter += weight;
            }
            return;
        }
        // ---------- NO PLANE ----------
        // else doesn't need else-brackets due to return statement
        if (cppf[1] > 0 || cppf[2] < 0) // "quadrant/octant conditions": is second point on right side of tangent plane?
            return;
        // compute "sppf" to check cylinder + sphere conditions
        const T S1 = std::sqrt(cppf[0] - cppf[1] * cppf[1]);
        const T S2 = std::sqrt(cppf[0] - cppf[2] * cppf[2]);
        const T S3 = std::sqrt(1. - cppf[3] * cppf[3]);
        const T cos_triangle = cppf[1] * cppf[2] * cppf[3] + S1 * S2 * cppf[3] + S1 * cppf[2] * S3 - cppf[1] * S2 * S3;
        const T cos_symmetry = S1 * S2 - cppf[1] * cppf[2];
        if (cos_symmetry >= cppf[0] || cos_triangle >= cppf[0]) // data badly conditioned
            return;
        // ---------- CONE ----------
        T dist = axis_dist(cppf);
//        Vec3 a = axis(cppf, pi, ni);
//        T sin_t = std::abs(a.dot(this->nr_));
//        dist = (1. - sin_t * sin_t) * dist / sin_t;
        if (!covt.not_valid(dist)) // only compute cone axis if distance is valid
            // no ClosenessWeights check necessary, since no equality constraints present for cone
            covt.vote(dist, axis(cppf, pi, ni), s.InterpolatedWeights);
        // ---------- CYLINDER OR SPHERE ----------
        if (cos_symmetry < this->cos_thresh2_ * cppf[0]) // normal projections have to be symmetric w.r.t distance vector
            return;
        // ---------- CYLINDER ----------
        if (!s.ClosenessWeights) {
            if (cos_symmetry > this->cos_thresh_ * cppf[0])
                cyvt.vote(radius(cppf), angle(ni), s.InterpolatedWeights);
        }
        else {
            const T bin_inv = this->angle_bin_inv_;
            T weight = (1. - .5 * std::acos(cos_symmetry / cppf[0]) * bin_inv);
            cyvt.vote(radius(cppf), angle(ni), weight, s.InterpolatedWeights);
        }
        // ---------- SPHERE ----------
        if (cos_triangle < this->cos_thresh_ * cppf[0]) // condition for sphere: normals and distance vector form triangle
            return;
        if (!s.ClosenessWeights) {
            if (cos_symmetry > this->cos_thresh_ * cppf[0] && cos_triangle > this->cos_thresh_half_ * cppf[0])
                svt.vote(radius(cppf), s.InterpolatedWeights);
        }
        else {
            const T bin_inv = this->angle_bin_inv_;
            T weight = (1. - .5 * std::acos(cos_symmetry / cppf[0]) * bin_inv);
            weight *= (1. - std::acos(cos_triangle / cppf[0]) * bin_inv);
            svt.vote(radius(cppf), weight, s.InterpolatedWeights);
        }
        return;
    }
    
    virtual void add_candidate(Vec3 pr, Vec3 nr) {
        using Vec2 = Eigen::Matrix<T, 2, 1>;
        VoteT max_count = this->min_votes_;
        int Obj = 0; // no candidate
        if (plane_counter > max_count) {
            max_count = plane_counter;
            Obj = 1; // plane
        }
        svt.set_threshold(max_count);
        std::pair<VoteT, T> el_sphere = svt.get_radius(this->settings_.MeanShift);
        if (el_sphere.first > max_count) {
            max_count = el_sphere.first;
            Obj = 2; // sphere
        }
        cyvt.set_threshold(max_count);
        std::pair<VoteT, Vec2> el_cyl = cyvt.get_radius_angle(this->settings_.MeanShift);
        if (el_cyl.first > max_count) {
            max_count = el_cyl.first;
            Obj = 3; // cylinder
        }
        covt.set_threshold(max_count);
        std::pair<VoteT, Vec4> el_cone = covt.get_dist_axis(this->settings_.MeanShift);
        if (el_cone.first > max_count) {
            max_count = el_cone.first;
            Obj = 4; // cone
        }
        switch (Obj) {
            case 0 :
                break;
            case 1 : {
                plane_candidates_.emplace(max_count, new Plane<T>(nr, -nr.dot(pr), pr));
                break;
            }
            case 2 : {
                T radius = el_sphere.second;
                Vec3 c = pr - radius * nr;
                sphere_candidates_.emplace(max_count, new Sphere<T>(c, radius, pr));
                break;
            }
            case 3 : {
                T radius = el_cyl.second[0];
                T angle  = el_cyl.second[1];
                Vec3 c = pr - radius * nr;
                Vec3 a = this->rot_.transpose() * Vec2(std::cos(angle), std::sin(angle));
                cyl_candidates_.emplace(max_count, new Cylinder<T>(c, a, radius, pr));
                break;
            }
            case 4 : {
                T dist = el_cone.second[0];
                Vec3 a(el_cone.second[1], el_cone.second[2], el_cone.second[3]);
                Vec3 c = pr + dist * (a / a.dot(nr) - nr);
                a = (a.dot(pr - c) > 0) ? a : -a;
                T sin_theta = -a.dot(nr);
//                T sin_theta = a.dot(nr);
//                a = sin_theta > 0 ? -a : a;
//                sin_theta = std::abs(sin_theta);
//                Vec3 c = pr - dist * (a + sin_theta * nr) / (1 - sin_theta * sin_theta);
                // discard if opening angle too small (i.e. close to cylinder), or too large (i.e. close to plane)
                // but should not actually be necessary since these options are already excluded --> TODO: investigate/debug!
                if (sin_theta > this->sin_thresh_ && sin_theta < this->cos_thresh2_)
                    cone_candidates_.emplace(max_count, new Cone<T>(c, a, sin_theta, pr));
                break;
            }
        }
    }

public:

    PrimitiveDetector(DetectorSettings settings = DetectorSettings(false, false, false, false, false, false), T Rmax = 1.) :
        PpfDetector<Pcd, T, VoteT>(settings),
        svt(SphereVotingTable<T, VoteT>(Rmax, this->dist_bin_)),
        cyvt(CylinderVotingTable<T, VoteT>(Rmax, this->dist_bin_, this->angle_bin_)),
        covt(ConeVotingTable<T, VoteT>(Rmax, this->dist_bin_, this->angle_bin_)) // Dmax = Rmax?!
    {
        this->settings_.PpfConditions = false; // not yet implemented
        this->settings_.ApproximatedConditions = false; // not yet implemented
        svt.set_threshold(this->min_votes_);
        cyvt.set_threshold(this->min_votes_);
        covt.set_threshold(this->min_votes_);
    }
    
    PrimitiveDetector(T diameter, bool is_scene, DetectorSettings settings = DetectorSettings(false, false, false, false, false, false)) :
        PpfDetector<Pcd, T, VoteT>(settings, diameter, is_scene),
        svt(SphereVotingTable<T, VoteT>(this->max_param_, this->dist_bin_)),
        cyvt(CylinderVotingTable<T, VoteT>(this->max_param_, this->dist_bin_, this->angle_bin_)),
        covt(ConeVotingTable<T, VoteT>(this->max_param_, this->dist_bin_, this->angle_bin_)) // Dmax = Rmax?!
    {
        this->settings_.PpfConditions = false; // not yet implemented
        this->settings_.ApproximatedConditions = false; // not yet implemented
        svt.set_threshold(this->min_votes_);
        cyvt.set_threshold(this->min_votes_);
        covt.set_threshold(this->min_votes_);
    }

    virtual ~PrimitiveDetector() {
        for (auto& el : plane_candidates_) {
            delete el.second;
        }
        for (auto& el : sphere_candidates_) {
            delete el.second;
        }
        for (auto& el : cyl_candidates_) {
            delete el.second;
        }
        for (auto& el : cone_candidates_) {
            delete el.second;
        }
    }
    
    virtual void cluster_candidates() {
        this->template cluster_candidates_base<Plane<T>>(plane_candidates_);
        this->template cluster_candidates_base<Sphere<T>>(sphere_candidates_);
        this->template cluster_candidates_base<Cylinder<T>>(cyl_candidates_);
        this->template cluster_candidates_base<Cone<T>>(cone_candidates_);
    }
    
    virtual void reset() {
        plane_counter = 0;
        svt.reset();
        cyvt.reset();
        covt.reset();
    }
    
    virtual std::vector<MyObject<T>*> candidate_vector(const T min_vote = 0, const T max_points = std::numeric_limits<T>::infinity()) const {
    
        std::multimap<T, MyObject<T>*, std::greater<T>> candidate_map;
        
        for (auto el : plane_candidates_) {
            candidate_map.emplace(std::sqrt(el.first), static_cast<MyObject<T>*>(el.second));
        }
        for (auto el : sphere_candidates_) {
            candidate_map.emplace(std::sqrt(el.first), static_cast<MyObject<T>*>(el.second));
        }
        for (auto el : cyl_candidates_) {
            candidate_map.emplace(std::sqrt(el.first), static_cast<MyObject<T>*>(el.second));
        }
        for (auto el : cone_candidates_) {
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
                  << sphere_candidates_.size() << " sphere candidates." << std::endl
                  << cyl_candidates_.size() << " cylinder candidates." << std::endl
                  << cone_candidates_.size() << " cone candidates." << std::endl;
    }

    virtual void write_candidates(std::ostream& os) const {
        os << "# planes:" << std::endl;
        for (const auto el : plane_candidates_) {
            os << el.first << "\t" << *(el.second) << std::endl;
        }
        os << "# spheres:" << std::endl;
        for (const auto el : sphere_candidates_) {
            os << el.first << "\t" << *(el.second) << std::endl;
        }
        os << "# cylinders:" << std::endl;
        for (const auto el : cyl_candidates_) {
            os << el.first << "\t" << *(el.second) << std::endl;
        }
        os << "# cones:" << std::endl;
        for (const auto el : cone_candidates_) {
            os << el.first << "\t" << *(el.second) << std::endl;
        }
    }
    
    virtual void write_results(const std::string folder, const std::string filename, const T min_vote = 0, const T max_points = std::numeric_limits<T>::infinity()) const {
    
        std::string prefix = filename.substr(0, filename.size() - 4);
        
        std::ofstream plane_file(folder + prefix + "_planes.txt");
        std::ofstream sphere_file(folder + prefix + "_spheres.txt");
        std::ofstream cyl_file(folder + prefix + "_cylinders.txt");
        std::ofstream cone_file(folder + prefix + "_cones.txt");
        
        std::vector<MyObject<T>*> candidates = candidate_vector(min_vote, max_points); // add second argument!!
        size_t instance = 1;
        for (const auto& obj : candidates) {        
            if (Plane<T>* p = dynamic_cast<Plane<T>*>(obj)) {
                plane_file << instance << "\t" << *p << std::endl;
            }
            else if (Sphere<T>* s = dynamic_cast<Sphere<T>*>(obj)) {
                sphere_file << instance << "\t" << *s << std::endl;
            }
            else if (Cylinder<T>* c = dynamic_cast<Cylinder<T>*>(obj)) {
                cyl_file << instance << "\t" << *c << std::endl;
            }
            else if (Cone<T>* c = dynamic_cast<Cone<T>*>(obj)) {
                cone_file << instance << "\t" << *c << std::endl;
            }
            ++instance;
        }
        
        plane_file.close();
        sphere_file.close();
        cyl_file.close();
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

#endif // PRIMITIVE_DETECTOR_H_
