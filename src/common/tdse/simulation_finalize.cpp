#include "common/tdse/simulation.h"

using namespace tdse;

void Simulation::WriteFinalState() const {
    _tdse_out->PushGroup("final_state");
    _tdse_out->WriteVector("wavefunction", _psi);
    _tdse_out->PopGroup();
}

void Simulation::Finish() {
    // ensure all resources are released
    _tdse_out = nullptr;
    _propagator = nullptr;
    _psi = nullptr;
    for (auto& p : _potentials)
        p = nullptr;
    for (auto& o : _observables)
        o = nullptr;
    for (auto& p : _pulses)
        p = nullptr;
}