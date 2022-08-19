#ifndef __ENERGY_RESOLVED_H__
#define __ENERGY_RESOLVED_H__

#include <string>
#include "common/ionization/ionization_task/ionization_task.h"
#include "common/utility/json.hpp"

class EnergyResolved : public IonizationTask::Register<EnergyResolved> {
    std::string _output_filename;
    double _emin, _emax, _estep;
public:
    static std::string GetName();
    static bool Validate(const nlohmann::json& energy_resolved);
    static IonizationTask::Ptr_t Create(const nlohmann::json& energy_resolved);

    void Execute();
};


#endif