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

#include <iostream>
#include <Eigen/Dense>
#include <ceres/ceres.h>
// folder includes
#include "TruncatedLoss.h"
// function includes
#include "pcd/PointCloudData.h"

template <class Pcd>
class ShapeOptimizer {

    const double median_mult_ = 4.;
    const double initial_scale_ = 0.05;

    // approximate median computation (enough for our purpose)
    static double median(std::vector<double>& vec) {
        size_t n = vec.size();
        if (n < 2) {
            return 0;
        }
        else {
            std::nth_element(vec.begin(), vec.begin() + n / 2, vec.end());
            return vec[n / 2];
        }
    }
    
    // compute right scale for ceres::LossFunction
    void recompute_scale() {
        std::vector<double> residuals, abs_residuals;
        ceres::Problem::EvaluateOptions options;
        options.apply_loss_function = false;
        problem_->Evaluate(options, nullptr, &residuals, nullptr, nullptr);

        if (residuals.size() < 3) {
            return;
        }

        for (const auto& r : residuals) {
            double r_abs = std::abs(r);
            if (r_abs < scale_)
                abs_residuals.push_back(r_abs);
        }

        scale_ = median_mult_ * median(abs_residuals);
    }    
    
protected:
    
    ceres::Problem* problem_;
    
    double scale_ = initial_scale_;
    
    virtual void add_residual(Eigen::Vector3f point) = 0;
    
    virtual void set_param() {}

public:

    bool optimize(const Pcd& pc) {
        
        // set up solver
        ceres::Solver::Options options;
        options.max_num_iterations = 5;
        options.linear_solver_type = ceres::DENSE_QR;
        options.minimizer_progress_to_stdout = false;//true;
        ceres::Solver::Summary summary;
          
        for (size_t outer=0; outer<25; ++outer) {
            
            // set up problem
            problem_ = new ceres::Problem;
            for (size_t idx=0; idx<pc.num_points(); ++idx) {
                add_residual(pc.point(idx));
            }
            set_param();
            
            // solve problem
            ceres::Solve(options, problem_, &summary);
//            std::cout << "Current scale: " << scale_ << std::endl;
//            std::cout << summary.BriefReport() << std::endl; // BriefReport or FullReport
            
            if (summary.num_successful_steps < 2) {
                delete problem_;
                break;
            }
            
            recompute_scale();
            
            delete problem_;
        }
        
        std::cout << "Scale: " << scale_ << ", convergence: " << (summary.termination_type == ceres::CONVERGENCE) << std::endl;
        
        if (scale_ > 0.04 || summary.termination_type == ceres::NO_CONVERGENCE) {
            scale_ = initial_scale_;
            return false;
        }

        scale_ = initial_scale_;
        return true;
    }

};
