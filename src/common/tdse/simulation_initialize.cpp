#include "common/tdse/simulation.h"
#include "common/utility/logger.h"
#include "common/utility/file_exists.h"
#include "common/utility/index_manip.h"
#include "common/file_io/io_factory.h"
#include "common/maths/math_factory.h"
#include "common/maths/math_algebra.h"
#include "common/bspline/bspline.h"
#include "common/system_state/system_state.h"

using namespace tdse;

void Simulation::ComputeDuration() {
    // find total length of simulation by finding the latest nonzero pulse.
    _tmax = 0; 
    for (auto& p : _pulses)
        _tmax = std::max(_tmax, p->delay + p->duration);
    _NT =  _tmax / GetTimestep() + 1;
}

void Simulation::InitializePsi() {
    _psi = maths::Factory::CreateVector(_dof);
    // if do_prop is set to false, we can assume the user
    // wants to just run the startup/shutdown observables.
    // Therefore load the most recent wave function if we can.
    if (!_do_propagate)       
        _restarting = true;
        
    // are we restarting?
    if (_restarting)
        _restarting = file_exists("TDSE.h5");

    // open checkpoint file
    _tdse_out = io::Factory::OpenHDF5("TDSE.h5", (_restarting ? 'a' : 'w'));

    if (_restarting) {          // if we are restarting make sure the input file matched
        if (!CompareTDSEH5wInput())
            exit(-1);           // failure
    } else 
        WriteParametersToTDSE();

    // if we are not restarting or there are no recent checkpoints
    if (!_restarting || !LoadLastCheckpoint()) {
        LoadInitialState();
        WriteInitialState();
    }
}

void Simulation::LoadInitialState() {
    LOG_INFO("Loading initial state...");
    _initialState.BuildInitialState(_psi);
    LOG_INFO("done...");
}
bool Simulation::LoadLastCheckpoint() {
    int it;
    LOG_INFO("Loading last checkpoint...");
    // what was the last successful checkpoint?
    _tdse_out->PushGroup("parameters");
    _tdse_out->ReadAttribute("last_checkpoint", &it);
    _tdse_out->PopGroup();

    if (it == -1) {             // ran but never had a checkpoint
        LOG_INFO("No checkpoints found...");
        return false;
    }
    _startIteration = it;
    // load that checkpoint
    _tdse_out->PushGroup("checkpoints");
    _tdse_out->ReadVector(std::to_string(_startIteration), _psi);
    _tdse_out->PopGroup();

    return true;
}


void Simulation::ComputeFields() {
    double dt = GetTimestep();
    for (auto& p : _pulses) {
        if (p->polarization_vector.x != 0. || p->minor_polarization_vector.x != 0.)
            _field[maths::X].resize(_NT);
        if (p->polarization_vector.y != 0. || p->minor_polarization_vector.y != 0.)
            _field[maths::Y].resize(_NT);
        if (p->polarization_vector.z != 0. || p->minor_polarization_vector.z != 0.)
            _field[maths::Z].resize(_NT);
    }

    if (_field[maths::X].size() > 0) {
        for (int it = 0; it < _NT; it++) {
            for (auto& p : _pulses)
                _field[maths::X][it] += p->A(it*dt).x;
        }
    }
    if (_field[maths::Y].size() > 0) {
        for (int it = 0; it < _NT; it++)
            for (auto& p : _pulses)
                _field[maths::Y][it] += p->A(it*dt).y;
    }
    if (_field[maths::Z].size() > 0) {
        for (int it = 0; it < _NT; it++)
            for (auto& p : _pulses)
                _field[maths::Z][it] += p->A(it*dt).z;
    }
}

void Simulation::WriteInitialState() const {
    _tdse_out->PushGroup("initial_state");
    _tdse_out->WriteVector("wavefunction", _psi);
    _tdse_out->PopGroup();
}

void Simulation::WriteParametersToTDSE() const {
    _tdse_out->PushGroup("parameters");
    _tdse_out->WriteAttribute("time_step", GetTimestep());
    _tdse_out->WriteAttribute("time_max", _tmax);
    _tdse_out->WriteAttribute("num_timesteps", _NT);
    _tdse_out->WriteAttribute("x_min", bspline::Basis::GetXmin());
    _tdse_out->WriteAttribute("x_max", bspline::Basis::GetXmax());
    _tdse_out->WriteAttribute("l_max", SystemState::GetBasisLmax());
    _tdse_out->WriteAttribute("m_max", SystemState::GetBasisMmax());
    _tdse_out->WriteAttribute("nodes", bspline::Basis::GetNumNodes());
    _tdse_out->WriteAttribute("order", bspline::Basis::GetOrder());
    _tdse_out->WriteAttribute("ecs_r0", bspline::Basis::GetECS_R0());
    _tdse_out->WriteAttribute("ecs_theta", bspline::Basis::GetECS_Theta());
    _tdse_out->WriteAttribute("last_checkpoint", -1);        // for later
    _tdse_out->PopGroup();
}

#define COMPARE_PARAM(x,y) if (x != y) { \
    LOG_CRITICAL("input "#x" does not match TDSE.h5. Expected " + std::to_string(x)); \
    return false; \
}

bool Simulation::CompareTDSEH5wInput() const {
    // read parameters and compare to input file
    double ecs_r0, ecs_theta;
    int order, nodes;
    int lmax, mmax;
    double xmax, xmin;
    double dt, tmax;
    int NT;

    _tdse_out->PushGroup("parameters");
    _tdse_out->ReadAttribute("time_step", &dt);
    _tdse_out->ReadAttribute("time_max", &tmax);
    _tdse_out->ReadAttribute("num_timesteps", &NT);
    _tdse_out->ReadAttribute("x_min", &xmin);
    _tdse_out->ReadAttribute("x_max", &xmax);
    _tdse_out->ReadAttribute("nodes", &nodes);
    _tdse_out->ReadAttribute("order", &order);
    _tdse_out->ReadAttribute("l_max", &lmax);
    _tdse_out->ReadAttribute("m_max", &mmax);
    _tdse_out->ReadAttribute("ecs_r0", &ecs_r0);
    _tdse_out->ReadAttribute("ecs_theta", &ecs_theta);
    _tdse_out->PopGroup();

    // make sure everything matches
    COMPARE_PARAM(dt,GetTimestep());
    COMPARE_PARAM(tmax,_tmax);
    COMPARE_PARAM(NT,_NT);
    COMPARE_PARAM(xmin,bspline::Basis::GetXmin());
    COMPARE_PARAM(xmax,bspline::Basis::GetXmax());
    COMPARE_PARAM(lmax,SystemState::GetBasisLmax());
    COMPARE_PARAM(mmax,SystemState::GetBasisMmax());
    COMPARE_PARAM(nodes,bspline::Basis::GetNumNodes());
    COMPARE_PARAM(order,bspline::Basis::GetOrder());
    COMPARE_PARAM(ecs_r0,bspline::Basis::GetECS_R0());
    COMPARE_PARAM(ecs_theta,bspline::Basis::GetECS_Theta());
    
    return true;
}
