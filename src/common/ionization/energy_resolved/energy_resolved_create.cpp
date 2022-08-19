#include "common/ionization/energy_resolved/energy_resolved.h"

IonizationTask::Ptr_t EnergyResolved::Create(const nlohmann::json& energy_resolved) {
    EnergyResolved* new_er = new EnergyResolved();

    auto& energy = energy_resolved["energy"];
    new_er->_output_filename = energy_resolved["filename"];
    new_er->_emin = energy["min"];
    new_er->_emax = energy["max"];
    new_er->_estep = energy["step"];

    return IonizationTask::Ptr_t(new_er);
}
