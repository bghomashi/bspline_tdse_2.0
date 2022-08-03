#ifndef __PETSC_EPS_H__
#define __PETSC_EPS_H__

#include "common/maths/eigen_solver.h"

class PetscEPS : public maths::IEigenSolver {
public:
    void Solve(const maths::Matrix A, const maths::Matrix S, int numVectors, double tol);
};

#endif