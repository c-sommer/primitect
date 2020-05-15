/**
GPL (v3+) License

This file is part of the code accompanying the paper
PrimiTect: Fast Continuous Hough Voting for Primitive Detection
by C. Sommer, Y. Sun, E. Bylow and D. Cremers,
accepted for publication in the IEEE International Conference on Robotics and Automation (ICRA) 2020.

Copyright (c) 2017, Christiane Sommer.
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

#ifndef TIMER_H_
#define TIMER_H_

// includes
#include <iostream>
#include <ctime>
#include <string>

class Timer {

private:
    clock_t start_time;
    clock_t end_time;
    double elapsed;

public:
    Timer() : start_time(0.), end_time(0.), elapsed(0.) {}
    ~Timer() {}
    
    void tic() {
        start_time = clock();
    }
    
    double toc(std::string s = "Time elapsed") {
        if (start_time!=0) {
            end_time = clock();
            elapsed = double(end_time-start_time) / CLOCKS_PER_SEC;
            print_time(s);
        }
        else
            std::cout << "Timer was not started, no time could be measured." << std::endl;
        return elapsed;
    }
    
    void print_time(std::string s = "Time elapsed") {        
        if (elapsed<1.)
            std::cout << "---------- " << s << ": " << 1000.*elapsed << "ms." << std::endl;
        else
            std::cout << "---------- " << s << ": " << elapsed << "s." << std::endl;
    }

};

#endif // TIMER_H_
