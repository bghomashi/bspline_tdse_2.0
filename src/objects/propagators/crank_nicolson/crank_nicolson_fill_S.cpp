
#include "objects/propagators/crank_nicolson/crank_nicolson.h"
#include "common/utility/index_manip.h"
#include "common/tdse/simulation.h"
#include "common/bspline/bspline.h"

using namespace std::complex_literals;
using namespace tdse;
using Maths = maths::Factory;
using namespace maths;

void CrankNicolson::FillOverlap(Matrix& S) {
    int N = bspline::Basis::GetNumBSplines();
    int order = bspline::Basis::GetOrder();
    int dof = Simulation::GetDOF();
    auto& Ms = Simulation::GetMs();
    auto& MRows = Simulation::GetMRows();

    std::vector<complex> overlapStore(N*N);
    for (int i = 0; i < N; i++) {
        for (int j = i; j < N; j++) {
            overlapStore[j + i*N] = overlapStore[i + j*N] = bspline::Basis::Integrate(i+1, j+1);
        }
    }
    S->FillBandedBlock(order-1, N, [=,&overlapStore](int row, int col) {
        int i, j, l1, l2, m1, m2;
        ILMFrom(row, i, l1, m1, N, Ms, MRows);
        ILMFrom(col, j, l2, m2, N, Ms, MRows);

        return overlapStore[i + j*N];
    });
}
