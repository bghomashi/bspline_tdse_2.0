#include "common/tise/tise.h"
#include "common/utility/logger.h"
#include "common/utility/to_lower.h"
#include "common/bspline/bspline.h"
#include "common/system_state/system_state.h"

using namespace tise;

 bool TISE::_Validate(const nlohmann::json& input) const {
    LOG_INFO("Validating input file.");

    if (!SystemState::Validate(input))
        return false;
    
    auto& eigen_state = input["eigen_state"];

    // --------- solver
    if (!eigen_state.contains("solver") || !eigen_state["solver"].is_string()) {
        LOG_CRITICAL("\"eigen_state\" must contain string entry: solver.");
        return false; 
    }
    if (ToLower(eigen_state["solver"]) != "slepc") {
        LOG_CRITICAL("eigen_state only supports \"SLEPc\" solver.");
        return false; 
    }
    
    // --------- tol
    if (eigen_state.contains("tol") && !eigen_state["tol"].is_number()) {
        LOG_CRITICAL("Optional entry \"tol\" in eigen_state must be a number.");
        return false; 
    }
    if (eigen_state.contains("tol") && eigen_state["tol"].get<double>() < 0) {
        LOG_CRITICAL("eigen_state \"tol\" must be positive.");
        return false; 
    }
    // --------- ecs
    if (eigen_state.contains("ecs_on") && !eigen_state["ecs_on"].is_boolean()) {
        LOG_CRITICAL("Optional entry \"ecs_on\" in eigen_state must be a boolean.");
        return false; 
    }

    // ---------- basis
    if (!bspline::Basis::Validate(input["basis"])) {
        LOG_CRITICAL("Failed to contruct \"basis\".");
        return false;
    }

    // ---------- continuum
    if (eigen_state.contains("continuum")) {
        auto& continuum = eigen_state["continuum"];
        if (continuum.contains("normalization") && !continuum["normalization"].is_string()) {
            LOG_CRITICAL("Optional entry \"eigen_state.continuum.normalization\" must be a string.");
            return false;
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
            if (!tdse::Potential::Validate(term["type"], term)) {
                LOG_CRITICAL("Failed to contruct \"potentials\".");
                return false;
            }
        }
    }

    return true;
}