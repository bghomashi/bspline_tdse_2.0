#include "common/ionization/ionization.h"
#include "common/utility/logger.h"
#include "common/system_state/system_state.h"
#include "common/bspline/bspline.h"

bool Ionization::Validate(const nlohmann::json& input) {
    LOG_INFO("Validating input file.");

    if (!SystemState::Validate(input))
        return false;
    
    if (!input.contains("ionization") || !input["ionization"].is_object()) {
        LOG_CRITICAL("Input file requires \"ionization\" object.");
        return false;
    }

    auto& ionization = input["ionization"];
    auto& eigen_state = input["eigen_state"];
    
    // --------- ecs
    if (eigen_state.contains("ecs_on") && !eigen_state["ecs_on"].is_boolean()) {
        LOG_CRITICAL("Optional entry \"ecs_on\" in eigen_state must be a boolean.");
        return false; 
    }
    // ---------------- basis ----------------
    if (!bspline::Basis::Validate(input["basis"])) {
        LOG_CRITICAL("Failed to contruct \"basis\".");
        return false;
    }

    // --------------- ionization observable ---------
    if (!input.contains("observables")) {
        LOG_CRITICAL("Input file requires \"observables\" object.");
        return false;
    }
    if (!input["observables"].contains("ionization")) {
        LOG_CRITICAL("\"Observables\" requires object entry: ionization.");
        return false;
    }
    if (!input["observables"]["ionization"].contains("filename") || !input["observables"]["ionization"]["filename"].is_string()) {
        LOG_CRITICAL("\"ionization\" observable requires string entry: filename.");
        return false;
    }

    // --------------- ionization tasks ----------
    for (auto& iot_pair : ionization.items()) {
        if (!IonizationTask::Validate(iot_pair.key(), iot_pair.value())) {
            LOG_CRITICAL("Failed to contruct \"observables\".");
            return false;
        }
    }
    
    return true;
}