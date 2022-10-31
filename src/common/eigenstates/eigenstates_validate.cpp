#include "common/eigenstates/eigenstates.h"
#include "common/utility/logger.h"
#include "common/system_state/system_state.h"
#include "common/bspline/bspline.h"

bool Eigenstates::Validate(const nlohmann::json& input) {
    LOG_INFO("Validating input file.");

    if (!SystemState::Validate(input))
        return false;
    
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
    
    return true;
}