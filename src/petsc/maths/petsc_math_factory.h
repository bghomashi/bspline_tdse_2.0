#ifndef __PETSC_MATH_FACTORY_H__
#define __PETSC_MATH_FACTORY_H__

#include "common/maths/math_factory.h"

class PetscMathFactory : public maths::IFactory {
    int _size, _rank;
public:
    bool Startup(int argc, char **args);
    void Shutdown();
    
    maths::Vector CreateVector(int N);
    void DestroyVector(maths::Vector& m);

    maths::Matrix CreateMatrix(int rows, int cols, int numBands);
    void DestroyMatrix(maths::Matrix& m);

    maths::GMRESSolver CreateGMRESSolver(int restart_iter = 500, int max_iter = 10000);
    void DestroyGMRESSolver(maths::GMRESSolver& m);

    maths::EigenSolver CreateEigenSolver();
    void DestroyEigenSolver(maths::EigenSolver& m);
};

#endif