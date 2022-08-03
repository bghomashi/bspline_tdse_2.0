#include "common/maths/eigen_solver.h"

using namespace maths;

IEigenSolver::IEigenSolver() : 
    _problemType(EigenProblemType::GNHEP) {}

const std::vector<complex>& IEigenSolver::GetEigenValues() const {
    return _values;
}
const std::vector<Vector>& IEigenSolver::GetEigenVectors() const {
    return _vectors;
}
void IEigenSolver::SetProblemType(EigenProblemType pt) {
    _problemType = pt;
}