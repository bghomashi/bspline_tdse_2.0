#include "common/tdse/simulation.h"
#include "common/utility/logger.h"
#include "common/bspline/bspline.h"
#include "common/system_state/system_state.h"

using namespace tdse;

void Simulation::_Load(const nlohmann::json& input) {
    SystemState::Load(input);

    _initialState.Load(input["initial_state"]);                 // parse input file (does not load actual state)
    _initialState._eigenStateFilename = SystemState::GetEigenStateFilename();
    
    bspline::Basis::Load(input["basis"]);

    if (input.contains("lasers") && input["lasers"].is_array()) {
        for (auto& pulse : input["lasers"])
            AddPulse(Pulse::Create(pulse["envelope"], pulse));
    }

    if (input.contains("observables")) {
        auto& observables = input["observables"];
        for (auto& obs_pair : observables.items()) {
            AddObservable(Observable::Create(obs_pair.key(), obs_pair.value()));
        }
    }

    if (input.contains("potentials")) {
        for (auto& pot : input["potentials"]) {
            AddPotential(Potential::Create(pot["type"], pot));
        }
    }
    
    CheckSymmetry();

    _checkpoints = input["propagator"]["checkpoint"].get<int>();
    if (input["propagator"].contains("restart"))
        _restarting = input["propagator"]["restart"];
    if (input["propagator"].contains("do_propagate"))
        _do_propagate = input["propagator"]["do_propagate"];

    _propagator = Propagator::Create(input["propagator"]["type"], input["propagator"]);
    _propagator->Initialize();
}