
#include "common/tise/tise.h"
#include "common/utility/banded_matrix.h"
#include "common/utility/logger.h"
#include "common/utility/profiler.h"
#include "common/utility/file_exists.h"
#include "common/file_io/io_factory.h"
#include "common/maths/math_algebra.h"
#include "common/maths/eigen_solver.h"
#include "common/bspline/bspline.h"
#include "common/system_state/system_state.h"

#include <string>
#include <sstream>
#include <complex>
#include <iomanip>

// in the axial potential case all the M's are decoupled. 
// The L is not a good quantum number (nbut M is). We expand in  
// sphericalharmonics anywa and have to diagonalize all the L's 
// together.The N*lmax x N*lmax matrices  will have 
// at most (2*order-1)*lmax bands.

using namespace tdse;
using namespace maths;

void tise::TISE::SolveAxial() {
    LOG_INFO("Solving the axial EVP...");
    
    int N = bspline::Basis::GetNumBSplines();
    int order = bspline::Basis::GetOrder();
    int lmax = SystemState::GetBasisLmax();
    int numBands = (2*order - 1)*(lmax+1);
    int numBoundNstates = SystemState::GetEigenStateBoundNmax();
    int numContinuumNstates = SystemState::GetEigenStateContinuumNmax();
    auto outputFilename = SystemState::GetEigenStateFilename();
    EigenProblemType pt = SystemState::GetProblemType();
    Matrix H0, S, temp;                                              // 3 * N*(lmax + 1) * (2*order-1)*(lmax + 1)*16   Bytes
    BandedMatrix kinBlockStore(N, 2*order-1), r2BlockStore(N, 2*order-1), overlap(N, 2*order-1);  // 3*N*(2*order - 1) Bytes

    // ---------------- create matrices ---------------
    LOG_INFO("Creating matrices...");
    H0 = maths::Factory::CreateMatrix(N*(lmax+1), N*(lmax+1), numBands);
    S = maths::Factory::CreateMatrix(N*(lmax+1), N*(lmax+1), numBands);
    temp = maths::Factory::CreateMatrix(N*(lmax+1), N*(lmax+1), numBands);
    
    // --------------- cache common kinetic energy factors 
    LOG_INFO("Caching kinetic energy...");
    for (int i = 0; i < N; i++) {
        for (int j = i; j < std::min(N, i+order); j++) {
            kinBlockStore(i, j) = kinBlockStore(j, i) = bspline::Basis::Integrate(i+1, j+1, 1,1) / 2.0;
            r2BlockStore(i, j) =  r2BlockStore(j, i) = bspline::Basis::Integrate(i+1, j+1, [] (complex r) {
                return 1./r/r;
            });
            overlap(i, j) =  overlap(j, i) = bspline::Basis::Integrate(i+1, j+1);
        }
    }
    // --------------- fill overlap matrix
    LOG_INFO("Building overlap matrix.");
    S->FillBandedBlock(order-1, N, [=](int row, int col) {
        int i, j;
        i = row % N;
        j = col % N;
        
        return overlap(i, j);
    });

    
    LOG_INFO("Building Hamiltonian...");
    LOG_INFO("inserting kinetic energy term");
    // --------------- fill all the L-blocks kinetic energies - BLOCK DIAGONAL
    H0->FillBandedBlock(order-1, N, [=](int row, int col) {
        int i, j, l;
        i = row % N;
        j = col % N;
        l = row / N;
        
        // kinetic energy and centrifugal term
        return kinBlockStore(i, j) + 0.5*l*(l+1.)*r2BlockStore(i, j);
    });


    // ------------- Add Potential Terms -----------------
    LOG_INFO("inserting potential term(s)");
    for (auto& p : _potentials) {
        temp->Zero();
        p->FillMatrix(temp, N);          // get potential matrix
        AXPY(H0, 1., temp);              // add potential to H0
    }

    // -------------- overwrite existing file ----------
    io::Factory::OpenHDF5(SystemState::GetEigenStateFilename(), 'w');
    
    // ------------ Solve for eigenstates --------------
    if (numBoundNstates > 0) {
        LOG_INFO("Solve for bound states");
        _eigensolver->SetProblemType(pt);
        ComputeAndOutputBoundStates(numBoundNstates, H0, S);
    }
    if (numContinuumNstates > 0) {
        LOG_INFO("Solve for continuum states");
        EigenProblemType pt = SystemState::GetProblemType();
        _eigensolver->SetProblemType(pt);
        ComputeAndOutputContinuumStates(numContinuumNstates, H0, S);
    }
}
