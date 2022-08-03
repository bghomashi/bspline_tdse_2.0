#ifndef __BSPLINE_ECS_H__
#define __BSPLINE_ECS_H__

// class:       bspline::ECS
// description: Convenience class used to encapsulate the ECS evaluation.
//              By setting r0, theta the grid can be computed with the methods.
//              
//              ** In the input file r0 is the choose ratio of the grid where 
//              the first ECS node is found.
//              i.e.    r0 = ECS_NODE/(xmax - xmin);
//              Setting r0 = 1  (the last node on the grid) and/or theta = 0 
//              turns the ECS off.

#include "common/maths/math_common.h"

namespace bspline {
    struct ECS {
        double r0;					// [0, 1]
        double theta;				// radians 

        maths::complex q(double x) const;
        maths::complex R(double x) const;
        double x(maths::complex R) const;
    };
}


#endif