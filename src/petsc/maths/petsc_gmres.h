#ifndef __PETSC_SOLVER_H__
#define __PETSC_SOLVER_H__

#include <petsc.h>
#include "common/maths/gmres_solver.h"

class PetscSolver : public maths::IGMRESSolver {
    KSPConvergedReason _reason;
    const char *_strreason;
public:
    KSP _petsc_ksp;          /* linear solver context */
    PC _petsc_pc;            /* preconditioner context */

    PetscSolver(int restart_iter = 500, int max_iter = 10000);
    ~PetscSolver();

    void SetBlockedPC(int blocks);
    bool Solve(const maths::Matrix A, const maths::Vector b, maths::Vector x);
};

#endif