#include "common/current_density/current_density.h"
#include "common/utility/logger.h"
#include "common/system_state/system_state.h"
#include "common/bspline/bspline.h"

bool CurrentDensity::Validate(const nlohmann::json& input) {
    LOG_INFO("Validating input file.");

    if (!SystemState::Validate(input))
        return false;
        
    return true;
}