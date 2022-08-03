
#include "objects/propagators/crank_nicolson/crank_nicolson.h"
#include "common/utility/banded_matrix.h"
#include "common/utility/index_manip.h"
#include "common/tdse/simulation.h"
#include "common/maths/math_algebra.h"
#include "common/bspline/bspline.h"
#include "common/system_state/system_state.h"

using namespace std::complex_literals;
using namespace tdse;
using Maths = maths::Factory;
using namespace maths;

void CrankNicolson::FillFieldFree(Matrix& H0) {
    int N = bspline::Basis::GetNumBSplines();
    int order = bspline::Basis::GetOrder();
    int mmax = SystemState::GetBasisMmax();
    int lmax = SystemState::GetBasisLmax();
    int dof = Simulation::GetDOF();
    auto& Ms = Simulation::GetMs();
    auto& MRows = Simulation::GetMRows();
    auto& potentials = Simulation::GetPotentials();
    int potSymmetry = Simulation::GetPotentialSymmetry();


    // TODO: decide on banded structure in non-central case
    // - for now assuming central
    Matrix temp = Maths::CreateMatrix(dof, dof, (potSymmetry == Symmetry::Central ? 1 : lmax+1)*(2*order-1));
    BandedMatrix kinBlockStore(N, 2*order-1), r2BlockStore(N, 2*order-1);

    // Fill H0 with kinetic energy
    // derivative part is always the same (only depends on i,j)
    for (int i = 0; i < N; i++) {
        for (int j = i; j < std::min(i+order, N); j++) {
            kinBlockStore(i, j) = kinBlockStore(j, i) = bspline::Basis::Integrate(i+1, j+1, 1,1) /2.;
            r2BlockStore (i, j) =  r2BlockStore(j, i) = bspline::Basis::Integrate(i+1, j+1, [] (complex r) {
                return 1./r/r;
            });
        }
    }
    H0->FillBandedBlock(order-1, N, [=](int row, int col) {
        int i, j, l1, l2, m1, m2;
        ILMFrom(row, i, l1, m1, N, Ms, MRows);
        ILMFrom(col, j, l2, m2, N, Ms, MRows);
        return kinBlockStore(i, j) + 0.5*l1*(l1+1.)*r2BlockStore(i, j);
    });

    for (auto& p : potentials) {
        p->FillMatrix(temp, N, Ms, MRows);  // get potential matrix
        AXPY(H0, 1., temp);                 // add potential to H0
    }
}