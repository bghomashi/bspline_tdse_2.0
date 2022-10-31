#ifndef __CURRENT_DENSITY_H__
#define __CURRENT_DENSITY_H__

#include <string>
#include "common/maths/math_common.h"
#include "common/utility/json.hpp"
#include "common/ionization/ionization_task/ionization_task.h"
#include "common/objects/pulse.h"


class CurrentDensity {
    std::string _current_filename;
public:
    double _time_step;
    int continuum_nmax;                     
    int lmax, mmax, dof;
    std::vector<Pulse::Ptr_t> _pulses;

    void AddPulse(Pulse::Ptr_t p);
    bool Validate(const nlohmann::json& input);
    void Load(const nlohmann::json& input);
    void PopulationsCentral();
    void Execute();
};




#endif
