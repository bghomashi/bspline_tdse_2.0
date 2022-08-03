#include "common/tdse/simulation.h"
#include "common/utility/logger.h"
#include "common/objects/pulse.h"
#include "common/utility/file_exists.h"
#include "common/system_state/system_state.h"
#include "common/bspline/bspline.h"

using namespace tdse;


bool Simulation::_Validate(const nlohmann::json& input) const {
    LOG_INFO("Validating input file.");
    
    if (!SystemState::Validate(input))
        return false;

    if (!bspline::Basis::Validate(input["basis"])) {
        LOG_CRITICAL("Failed to contruct \"basis\".");
        return false;
    }

    // ----------- has initial state -----------
    if (!input.contains("initial_state") || !input["initial_state"].is_array()) {
        LOG_CRITICAL("input file must contain array entry: initial_state");
        return false;
    }
    if (!InitialState::Validate(input["initial_state"])) {
        LOG_CRITICAL("Failed to contruct \"initial_state\".");
        return false;
    }

    // ----------- propagator -----------
    if (!(input.contains("propagator") && input["propagator"].is_object())) {
        LOG_CRITICAL("input file must contain object entry: propagator");
        return false;
    } else {
        auto propagator = input["propagator"];
        if (!propagator.contains("type") || !propagator["type"].is_string()) {
            LOG_CRITICAL("propagator must contain string entry: type");
            return false;
        }

        if (!(propagator.contains("checkpoint") && propagator["checkpoint"].is_number())) {
            LOG_CRITICAL("propagator must contain number entry: checkpoint");
            return false;
        }

        if (propagator.contains("restart") && !propagator["restart"].is_boolean()) {
            LOG_CRITICAL("optional entry \"restart\" in propagator must be a boolean.");
            return false;
        }
        if (propagator.contains("do_propagate") && !propagator["do_propagate"].is_boolean()) {
            LOG_CRITICAL("optional entry \"do_propagate\" in propagator must be a boolean.");
            return false;
        }
        if (!Propagator::Validate(propagator["type"], propagator)) {
            LOG_CRITICAL("Failed to contruct \"propagator\".");
            return false;
        }
    }
    
    // ----------- lasers are optional -----------
    if (input.contains("lasers") && input["lasers"].is_array()) {
        auto& lasers = input["lasers"];

        for (auto& pulse : lasers) {
            if (!pulse.contains("envelope") || !pulse["envelope"].is_string()) {
                LOG_CRITICAL("each \"laser\" must contain string entry: envelope");
                return false;
            }

            if (!Pulse::Validate(pulse["envelope"], pulse)) {
                LOG_CRITICAL("Failed to contruct \"pulse\".");
                return false;    
            }
        }
    }

    // ----------- observables are optional -----------
    if (input.contains("observables")) {
        if (!input["observables"].is_object()) {
            LOG_CRITICAL("Optional entry \"observables\" must be an object.");
            return false;
        }
        auto& observables = input["observables"];
        for (auto& obs_pair : observables.items()) {
            if (!Observable::Validate(obs_pair.key(), obs_pair.value())) {
                LOG_CRITICAL("Failed to contruct \"observables\".");
                return false;
            }
        }
    }

    // ----------- potentials are optional -----------
    if (input.contains("potentials")) {
        if (!input["potentials"].is_array()) {
            LOG_CRITICAL("Optional entry \"potentials\" must be an array.");
            return false;
        }
        auto& potentials = input["potentials"];
        for (auto& term : potentials) {
            if (!(term.contains("type") && term["type"].is_string())) {
                LOG_CRITICAL("\"potential\" must contain string entry: type.");
                return false;
            }
            if (!Potential::Validate(term["type"], term)) {
                LOG_CRITICAL("Failed to contruct \"potentials\".");
                return false;
            }
        }
    }
            
    return true;
}