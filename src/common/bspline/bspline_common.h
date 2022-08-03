#ifndef __BSPLINE_COMMON_H__
#define __BSPLINE_COMMON_H__

// This file define a few constants used when generating the bspline basis.

#include "common/maths/math_common.h"

namespace bspline {
    constexpr double NOD_THRESHOLD = 1.E-15;								// grid spacing that is considered zero
    constexpr maths::complex I = maths::complex(0,1);
    constexpr double GAUSS_POINTS_CONVERGENCE_THRESHOLD = 2.E-15;	
    constexpr int NGAUSS = 64;	
    constexpr int DIMFACT = 127;	
}

#endif