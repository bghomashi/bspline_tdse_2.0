#include "common/tdse/simulation.h"
#include "common/system_state/system_state.h"
#include "common/utility/logger.h"
#include "common/bspline/bspline.h"

using namespace tdse;

void Simulation::CheckSymmetry() {
    int mmax = SystemState::GetBasisMmax();
    int lmax = SystemState::GetBasisLmax();
    int N = bspline::Basis::GetNumBSplines();
    // With the laser there is at most axial symmetry
    // ------------- potential symmetries ------------
    _potentialSymmetry = Symmetry::Central;
    _laserAxial = true;
    for (const auto& p : _potentials) {
        if (!p->isCentral())                    // l-blocks are not coupled
            _potentialSymmetry &= ~(1UL << 1);  // unset central bit
        if (!p->isAxial())                      // l-blocks are all coupled
            _potentialSymmetry &= ~(1UL << 0);  // unset axial bit   
    }

    // ------------- laser symmetries ------------
    for (const auto& p : _pulses) {
        _pol[maths::X] = p->polarization_vector.x || p->minor_polarization_vector.x;
        _pol[maths::Y] = p->polarization_vector.y || p->minor_polarization_vector.y;
        _pol[maths::Z] = p->polarization_vector.z || p->minor_polarization_vector.z;
    }
    if (_pol[maths::X] || _pol[maths::Y])
        _laserAxial = false;

    if (_laserAxial) {       // the Ms only come from the initial state and are not coupled
        for (auto& state : _initialState._initialState)
            _Ms.push_back(state.m);
    
        // sort and remove duplicates
        std::sort(_Ms.begin(), _Ms.end());
        _Ms.erase( std::unique(_Ms.begin(), _Ms.end() ), _Ms.end() );

        // still each m-block is independent so each row of interaction matrix
        // has 6*order - 3 elements
    } else {                // no symmetry - either due to potential or laser
        // we need all the m-values anyway
        _Ms.reserve(2*mmax+1);
        for (int m = -mmax; m <= mmax; m++)
            _Ms.push_back(m);

        // every m-block is coupled to the one before and after it
        // every row of the interaction has 12*_order - 6 elements
    }

    
    // ------------- tabulate m-block rows ---------------
    // count degrees of freedom and 
    // generate row table for m's
    _mRows.reserve(_Ms.size());
    _dof = 0;
    for (const int m : _Ms) {
        _mRows.push_back(_dof);
        if (_potentialSymmetry == Symmetry::Central)
            _dof += N*(lmax - std::abs(m) + 1);
        else if (_potentialSymmetry == Symmetry::Axial)
            _dof += N*(lmax + 1);
    }
}