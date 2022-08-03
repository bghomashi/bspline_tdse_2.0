#include "common/tdse/simulation.h"
#include "common/utility/logger.h"
#include "common/utility/profiler.h"


using namespace tdse;

void Simulation::_Execute() {
    double t = 0;
    int start_iteration = 0;

    // find total length of simulation by finding the latest nonzero pulse.
    LOG_INFO("Calculating cummulative field...");
    ComputeDuration();
    // calculate the field for each timestep in advance
    ComputeFields();
    LOG_INFO("Initializing wavefunction...");
    InitializePsi();

    // at this point TDSE.h5 is open.
    

    LOG_INFO("Executing Observable::Startup() ...");
    // allow observables to initialize
    // - always run start
    for (auto& obs : _observables)
        obs->Startup(start_iteration);

    if (!_do_propagate) {
        LOG_INFO("Not propagation...");
        _tdse_out = nullptr;

        // allow observables to complete
        LOG_INFO("Executing Observable::Shutdown() ...");
        for (auto& obs : _observables)
            obs->Shutdown();

        _propagator->Finish();          // make sure to cleanup anything that was initialized

        // close TDSE.h5
        _tdse_out = nullptr;

        return;
    }

    // OTHERWISE

    LOG_INFO("Beginning propagation...");
    // - propagate
    // - compute observables at the requested timesteps
    // - check points
    DoTimeSteps();

    WriteFinalState();

    _tdse_out = nullptr;

    // allow observables to complete
    LOG_INFO("Executing Observable::Shutdown() ...");
    for (auto& obs : _observables)
        obs->Shutdown();

    _propagator->Finish();          // make sure to cleanup anything that was initialized

    Finish();
}

void Simulation::DoTimeSteps() {
    LOG_INFO("Propagation for " + std::to_string(_NT) + " timesteps...");
    Profile::Push("Total time stepping");

    for (int it = _startIteration; it < _NT; it++) {
        if (!_propagator->DoStep(it))    // propagate
            break;        
        
        DoCheckpoint(it);               // checkpoint
        DoObservables(it);              // give observables a chance to execute
    }
    Profile::Pop("Total time stepping");
}

void Simulation::DoCheckpoint(int it) {
    if ((_checkpoints != 0) && (it % _checkpoints == 0)) {
        LOG_INFO("iteration: " + std::to_string(it) + "/" + std::to_string(_NT));
        for (auto& obs : _observables)
            obs->Flush();                           // dumb all the observables

        // dump psi to file
        _tdse_out->PushGroup("checkpoints");
        _tdse_out->WriteVector(std::to_string(it), _psi);
        _tdse_out->PopGroup();

        // write last checkpoint iteration
        _tdse_out->PushGroup("parameters");
        _tdse_out->WriteAttribute("last_checkpoint", it);           // should update
        _tdse_out->PopGroup();
    }
}

void Simulation::DoObservables(int it) {
    for (auto& obs : _observables)
        obs->DoObservable(it);
}