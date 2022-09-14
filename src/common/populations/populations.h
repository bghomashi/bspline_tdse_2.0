#ifndef __IONIZATION_H__
#define __IONIZATION_H__

#include <string>
#include "common/maths/math_common.h"
#include "common/utility/json.hpp"
#include "common/ionization/ionization_task/ionization_task.h"


class Populations {
    bool _ecs_on;
    std::string _populations_filename, _eigenstate_filename;
public:
    int continuum_nmax;                     
    int lmax, mmax, dof;

    bool Validate(const nlohmann::json& input);
    void Load(const nlohmann::json& input);
    void PopulationsCentral();
    void Execute();
};




#endif
