#ifndef __EIGEN_SOLVER_H__
#define __EIGEN_SOLVER_H__

#include "common/maths/math_common.h"
#include <vector>

namespace maths {
    enum struct EigenProblemType {
        GHEP,
        GNHEP,
        PGNHEP
    };
    class IEigenSolver {
    protected:
        std::vector<complex> _values;
        std::vector<Vector> _vectors;
        EigenProblemType _problemType;
    public:


        IEigenSolver();

        virtual void Solve(const Matrix A, const Matrix S, int numVectors, double tol) = 0;

        const std::vector<complex>& GetEigenValues() const;
        const std::vector<Vector>& GetEigenVectors() const;

        void SetProblemType(EigenProblemType pt);
    };
}

#endif