#ifndef __IONIZATION_H__
#define __IONIZATION_H__

#include <string>
#include "common/maths/math_common.h"
#include "common/utility/json.hpp"


class Eigenstates {
    bool _ecs_on;
    std::string _eigenstate_filename;
public:

    bool Validate(const nlohmann::json& input);
    void Load(const nlohmann::json& input);
    void PopulationsCentral();
    void Execute();
};




#endif
