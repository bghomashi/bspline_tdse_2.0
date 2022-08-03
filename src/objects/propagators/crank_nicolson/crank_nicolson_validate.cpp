#include "objects/propagators/crank_nicolson/crank_nicolson.h"
#include "common/utility/logger.h"

std::string CrankNicolson::GetName() {
    return "crank_nicolson";
}

tdse::Propagator::Ptr_t CrankNicolson::Create(const nlohmann::json& propagator) {
    auto prop = new CrankNicolson();

    prop->_dt = propagator["time_step"];

    return tdse::Propagator::Ptr_t(prop);
}
bool CrankNicolson::Validate(const nlohmann::json& propagator) {
    if (!(propagator.contains("time_step") && propagator["time_step"].is_number())) {
        LOG_CRITICAL("propagator must contain number entry: time_step");
        return false;
    }
    return true;
}