#include "common/tdse/simulation.h"
#include "common/utility/logger.h"
#include "common/utility/to_lower.h"

using namespace tdse;

bool Simulation::ValidateBasis(const nlohmann::json& basis) const {
    if (!(basis.contains("lmax") && basis["lmax"].is_number())) {
        LOG_CRITICAL("basis must contain number entry: lmax");
        return false;
    }
    if (!(basis.contains("mmax") && basis["mmax"].is_number())) {
        LOG_CRITICAL("basis must contain number entry: mmax");
        return false;
    }
    if (basis["mmax"] > basis["lmax"]) {
        LOG_CRITICAL("mmax must be less than or equal to lmax.");
        return false;
    }
    return true;
}