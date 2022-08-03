#ifndef __GMRES_SOLVER_H__
#define __GMRES_SOLVER_H__

#include "common/maths/math_common.h"

namespace maths {
    class IGMRESSolver {
    public:
        virtual bool Solve(const Matrix A, const Vector b, Vector x) = 0;
        virtual void SetBlockedPC(int blocks) = 0;
    };
}


#endif