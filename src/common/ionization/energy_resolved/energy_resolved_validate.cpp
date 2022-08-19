#include "common/ionization/energy_resolved/energy_resolved.h"
#include "common/utility/logger.h"

bool EnergyResolved::Validate(const nlohmann::json& energy_resolved) {
    if (!energy_resolved.contains("filename") || !energy_resolved["filename"].is_string()) {
        LOG_CRITICAL("\"energy-resolved\" must contain string entry: filename");
        return false;
    }
    if (!energy_resolved.contains("energy") || !energy_resolved["energy"].is_object()) {
        LOG_CRITICAL("\"energy-resolved\" must contain object entry: energy");
        return false;
    }
    auto& energy = energy_resolved["energy"];
    if (!energy.contains("min") || !energy["min"].is_number()) {
        LOG_CRITICAL("\"energy\" must contain number entry: min");
        return false;
    }
    if (!energy.contains("max") || !energy["max"].is_number()) {
        LOG_CRITICAL("\"energy\" must contain number entry: max");
        return false;
    }
    if (!energy.contains("step") || !energy["step"].is_number()) {
        LOG_CRITICAL("\"energy\" must contain number entry: step");
        return false;
    }
    return true;
}

std::string EnergyResolved::GetName() {
    return "energy_resolved";
}