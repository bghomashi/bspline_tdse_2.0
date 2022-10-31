
#include "objects/propagators/crank_nicolson/crank_nicolson.h"
#include <limits>

#include <iostream>
#include <sstream>
#include <string>
#include <algorithm>

#include "common/utility/index_manip.h"
#include "common/utility/logger.h"
#include "common/utility/profiler.h"
#include "common/maths/math_factory.h"
#include "common/maths/math_algebra.h"
#include "common/tdse/simulation.h"
#include "common/bspline/bspline.h"
#include "common/system_state/system_state.h"

#include "petsc/maths/petsc_common.h"

using namespace std::complex_literals;
using namespace tdse;
using Maths = maths::Factory;
using namespace maths;

CrankNicolson::CrankNicolson() {
}
void CrankNicolson::Initialize() {
    ProfilerPush();

    int N = bspline::Basis::GetNumBSplines();
    int order = bspline::Basis::GetOrder();
    int mmax = SystemState::GetBasisMmax();
    int lmax = SystemState::GetBasisLmax();

    int dof = Simulation::GetDOF();
    double dt = Simulation::GetTimestep();
    auto pol = Simulation::GetPolarization();
    int potSymmetry = Simulation::GetPotentialSymmetry();
    bool laserSymmetry = Simulation::GetLaserSymmetry();

    int bandsPerBlocks = 2*order - 1;
    int HIxBandsPerRow = 0;
    int HIyBandsPerRow = 0;
    int HIzBandsPerRow = 0;
    int U0BandsPerRow = 0;
    int H0BandsPerRow  = (potSymmetry == Symmetry::Central ? 1 : lmax+1) * bandsPerBlocks;

    if (pol[maths::X])
        HIxBandsPerRow = 4*bandsPerBlocks;
    if (pol[maths::Y])
        HIyBandsPerRow = 4*bandsPerBlocks;
    if (pol[maths::Z])
        HIzBandsPerRow = 2*bandsPerBlocks;

    U0BandsPerRow += H0BandsPerRow + HIxBandsPerRow + HIyBandsPerRow + HIzBandsPerRow;

    double memory = 0;
    memory += H0BandsPerRow*dof;            // memory for H0 matrix
    memory += bandsPerBlocks*dof;           // memory for S matrix (block diagonal)
    memory += HIxBandsPerRow*dof;           // memory for HIx matrix
    memory += HIyBandsPerRow*dof;           // memory for HIy matrix
    memory += HIzBandsPerRow*dof;           // memory for HIz matrix
    memory += 4*U0BandsPerRow*dof;          // memory for 4*U0 matrix
    memory += 2*dof;                        // wavefunctions
    memory = memory*16./1024./1024./1024.;  // in gigabytes

    LOG_INFO("Initialize TDSE");
    // ---------------------------------------------------------------
    // initialize field free, overlap, and interaction matrices
    LOG_INFO("Estimated memory required: " + std::to_string(memory) + " GB.");
    Log::Flush();

    LOG_INFO("Allocating space...");
    Matrix H0 = Maths::CreateMatrix(dof, dof, H0BandsPerRow);
    Matrix S  = Maths::CreateMatrix(dof, dof, bandsPerBlocks);
    
    if (pol[X])
        _HI[X] = Maths::CreateMatrix(dof, dof, HIxBandsPerRow);
    if (pol[Y])
        _HI[Y] = Maths::CreateMatrix(dof, dof, HIyBandsPerRow);
    if (pol[Z])
        _HI[Z] = Maths::CreateMatrix(dof, dof, HIzBandsPerRow);
    
    _U0p = Maths::CreateMatrix(dof, dof, U0BandsPerRow);
    _U0m = Maths::CreateMatrix(dof, dof, U0BandsPerRow);
    _Up = Maths::CreateMatrix(dof, dof, U0BandsPerRow);
    _Um = Maths::CreateMatrix(dof, dof, U0BandsPerRow);
    _psi_temp = Maths::CreateVector(dof);        // storage used to hold intermediate psi during propagation

    LOG_INFO("Filling intermediate matrices...");
    FillFieldFree(H0);
    FillOverlap(S);

    if (pol[X])
        FillInteractionX(_HI[X]);
    if (pol[Y])
        FillInteractionY(_HI[Y]);
    if (pol[Z])
        FillInteractionZ(_HI[Z]);

    // ---------------------------------------------------------------
    // Initialize the static propagator matrices (U0+/-)
    LOG_INFO("Building propagator matrix...");

    FillU0(_U0p);
    FillU0(_U0m);

    _U0p->AssembleBegin();
    _U0p->AssembleEnd();
    _U0m->AssembleBegin();
    _U0m->AssembleEnd();


    // add overlap (and zero off-diagonal blocks)
    AYPX(_U0p, 0, S);
    AYPX(_U0m, 0, S);
    // add fieldfree hamiltonian
    AXPY(_U0p, 0.5i*dt, H0);
    AXPY(_U0m, -0.5i*dt, H0);

    _Up->Duplicate(_U0p);
    _Um->Duplicate(_U0m);
    
    //-----------------------------------------------
    // Create solver
    _solver = Maths::CreateGMRESSolver();
    _solver->SetBlockedPC(N);
    //-----------------------------------------------

    LOG_INFO("Crank-Nicolson initialization complete.");

    ProfilerPop();

}

void CrankNicolson::Finish() {
    _U0p = nullptr;
    _U0m = nullptr;
    _HI[X] = nullptr;
    _HI[Y] = nullptr;
    _HI[Z] = nullptr;
    _Up = nullptr;
    _Um = nullptr;
    _psi_temp = nullptr;
    _solver = nullptr;
}
bool CrankNicolson::DoStep(int it) {
    auto& psi = Simulation::GetPsi();
    auto fields = Simulation::GetFields();
    
    _t = _dt*it;

    _Up->Copy(_U0p);
    _Um->Copy(_U0m);

    for (int xn = X; xn <= Z; xn++) {
        if (_HI[xn]) {
            // AXPY(_Up, 0.5i*_dt*(-1.i), _HI[xn]);
            // AXPY(_Um, -0.5i*_dt*(-1.i), _HI[xn]);
            AXPY(_Up, 0.5i*_dt*(-1.i*fields[xn][it]), _HI[xn]);
            AXPY(_Um, -0.5i*_dt*(-1.i*fields[xn][it]), _HI[xn]);
        }
    }
    
    _Um->AssembleBegin();
    _Um->AssembleEnd();
    // MatView(std::dynamic_pointer_cast<PetscMatrix>(_Um)->_petsc_mat, 0); exit(0);


    Mult(_Um, psi, _psi_temp);
    if (!_solver->Solve(_Up, _psi_temp, psi)) {
        std::cout << "divergence!" << std::endl;
        return false;               // failure
    }

    return true;
}

