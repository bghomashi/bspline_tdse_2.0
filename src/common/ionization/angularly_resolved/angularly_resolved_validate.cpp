#include "common/ionization/angularly_resolved/angularly_resolved.h"
#include "common/utility/logger.h"

bool AngularlyResolved::Validate(const nlohmann::json& angularly_resolved) {
    if (!angularly_resolved.contains("filename") || !angularly_resolved["filename"].is_string()) {
        LOG_CRITICAL("\"ionization.angularly_resolved\" must contain string entry: filename");
        return false;
    }
    if (!angularly_resolved.contains("energy") || !angularly_resolved["energy"].is_object()) {
        LOG_CRITICAL("\"ionization.angularly_resolved\" must contain object entry: energy");
        return false;
    }
    if (!angularly_resolved.contains("x") || !angularly_resolved["x"].is_object()) {
        LOG_CRITICAL("\"ionization.angularly_resolved\" must contain object entry: x");
        return false;
    }
    if (!angularly_resolved.contains("z") || !angularly_resolved["z"].is_object()) {
        LOG_CRITICAL("\"ionization.angularly_resolved\" must contain object entry: z");
        return false;
    }
    auto& energy = angularly_resolved["energy"];
    if (!energy.contains("min") || !energy["min"].is_number()) {
        LOG_CRITICAL("\"ionization.angularly_resolved.energy\" must contain number entry: min");
        return false;
    }
    if (!energy.contains("max") || !energy["max"].is_number()) {
        LOG_CRITICAL("\"ionization.angularly_resolved.energy\" must contain number entry: max");
        return false;
    }
    if (!energy.contains("step") || !energy["step"].is_number()) {
        LOG_CRITICAL("\"ionization.angularly_resolved.energy\" must contain number entry: step");
        return false;
    }
    auto& x = angularly_resolved["x"];
    if (!x.contains("min") || !x["min"].is_number()) {
        LOG_CRITICAL("\"ionization.angularly_resolved.x\" must contain number entry: min");
        return false;
    }
    if (!x.contains("max") || !x["max"].is_number()) {
        LOG_CRITICAL("\"ionization.angularly_resolved.x\" must contain number entry: max");
        return false;
    }
    if (!x.contains("step") || !x["step"].is_number()) {
        LOG_CRITICAL("\"ionization.angularly_resolved.x\" must contain number entry: step");
        return false;
    }
    auto& z = angularly_resolved["z"];
    if (!z.contains("min") || !z["min"].is_number()) {
        LOG_CRITICAL("\"ionization.angularly_resolved.z\" must contain number entry: min");
        return false;
    }
    if (!z.contains("max") || !z["max"].is_number()) {
        LOG_CRITICAL("\"ionization.angularly_resolved.z\" must contain number entry: max");
        return false;
    }
    if (!z.contains("step") || !z["step"].is_number()) {
        LOG_CRITICAL("\"ionization.angularly_resolved.z\" must contain number entry: step");
        return false;
    }
    return true;
}

std::string AngularlyResolved::GetName() {
    return "angularly_resolved";
}