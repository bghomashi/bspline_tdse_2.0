#ifndef __GSL_FIT_SIN_H__
#define __GSL_FIT_SIN_H__


#include <vector>


int fit_sin(
    const std::vector<double>& initial_guess,
    const std::vector<double>& ts,
    const std::vector<double>& ys,
    const std::vector<double>& constants,
    std::vector<double>& out);


#endif