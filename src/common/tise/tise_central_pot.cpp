
#include "common/tise/tise.h"
#include "common/utility/banded_matrix.h"
#include "common/utility/logger.h"
#include "common/utility/profiler.h"
#include "common/utility/file_exists.h"
#include "common/tdse/simulation.h"
#include "common/file_io/io_factory.h"
#include "common/maths/math_algebra.h"
#include "common/maths/eigen_solver.h"
#include "common/bspline/bspline.h"
#include "common/system_state/system_state.h"

#include <string>
#include <sstream>
#include <complex>
#include <iomanip>

// in the central potential case all the L's and M's are
// decoupled. We can diagonalize each L-block separately.
// The matrices are will have 2*order-1 bands. 


using namespace tdse;
using namespace maths;

void tise::TISE::SolveCentral() {
    int N = bspline::Basis::GetNumBSplines();
    int order = bspline::Basis::GetOrder();
    int numBands = 2*order - 1;
    int numBoundNstates = SystemState::GetEigenStateBoundNmax();
    int numBoundLstates = SystemState::GetEigenStateBoundLmax();
    int numContinuumNstates = SystemState::GetEigenStateContinuumNmax();
    int numContinuumLstates = SystemState::GetEigenStateContinuumLmax();
    int maxLstates = std::max(numBoundLstates, numContinuumLstates);
    BandedMatrix kinBlockStore(N, 2*order-1), r2BlockStore(N, 2*order-1);       // 2*N*(2*order-1)*16 Bytes
    Matrix H0, S;                                                               // 2*N*(2*order-1)*16 Bytes
    std::vector<Matrix> pot(_potentials.size());                                // _potential.size()*N*(2*order-1)*16 Bytes
    
    // ---------------- create matrices ---------------
    LOG_INFO("Allocating matrices...");
    H0 = maths::Factory::CreateMatrix(N, N, numBands);
    S = maths::Factory::CreateMatrix(N, N, numBands);
    for (auto& p : pot)
        p = maths::Factory::CreateMatrix(N, N, numBands);
    
    // --------------- cache common kinetic energy factors 
    LOG_INFO("Caching kinetic energy...");
    for (int i = 0; i < N; i++) {
        for (int j = i; j < std::min(N, i+order); j++) {
            kinBlockStore(i, j) = kinBlockStore(j, i) = bspline::Basis::Integrate(i+1, j+1, 1,1) / 2.0;
            r2BlockStore(i, j) =  r2BlockStore(j, i) = bspline::Basis::Integrate(i+1, j+1, [] (complex r) {
                return 1./r/r;
            });
        }
    }
 
    // --------------- fill overlap matrix
    LOG_INFO("Building overlap matrix.");
    S->FillBandedBlock(order-1, N, [=](int row, int col) {
        int i, j;
        i = row % N;
        j = col % N;
        
        return bspline::Basis::Integrate(i+1, j+1);
    });

    // -------------- overwrite existing file ----------
    io::Factory::OpenHDF5(SystemState::GetEigenStateFilename(), 'w');
    
    // -------------- fill potential matrices -------------
    for (int i = 0; i < _potentials.size(); i++)
        _potentials[i]->FillMatrix(pot[i], N);

    for (int l = 0; l <= maxLstates; l++) {
        LOG_INFO("Building Hamiltonian matrix (l=" + std::to_string(l) + ").");
        // ---------------- Fill H0 with kinetic energy
        H0->FillBandedBlock(order-1, N, [&](int row, int col) {
            int i, j;
            i = row % N;
            j = col % N;
            
            // kinetic energy and centrifugal term
            return kinBlockStore(i, j) + 0.5*l*(l+1.)*r2BlockStore(i, j);
        });

        // ------------- Add Potential Terms -----------------
        for (int i = 0; i < _potentials.size(); i++)
            AXPY(H0, 1., pot[i]);

        // ---------------- compute bound states --------------
        if (numBoundNstates > 0 && l <= numBoundLstates) {
            // ------------ Solve for eigenstates --------------
            LOG_INFO("Solve for bound states...");
            EigenProblemType pt = SystemState::GetProblemType();
            _eigensolver->SetProblemType(pt);
            ComputeAndOutputBoundStates(numBoundNstates - l, H0, S, l);
        }
        // ---------------- compute continuum states --------------
        if (numContinuumNstates > 0 && l <= numContinuumLstates) {
            LOG_INFO("Solve for continuum states...");
            EigenProblemType pt = SystemState::GetProblemType();
            _eigensolver->SetProblemType(pt);
            ComputeAndOutputContinuumStates(numContinuumNstates, H0, S, l);
        }
    }
}