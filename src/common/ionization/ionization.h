#ifndef __IONIZATION_H__
#define __IONIZATION_H__

#include <string>
#include "common/maths/math_common.h"
#include "common/utility/json.hpp"
#include "common/ionization/ionization_task/ionization_task.h"

struct Population {
    double energy;
    double population;
    double phaseshift;
    maths::complex amplitude;
};

class Ionization {
    bool _ecs_on;
    std::string _ionization_filename, _eigenstate_filename;
    std::vector<IonizationTask::Ptr_t> _ionization_tasks;


public:
    int continuum_nmax;                     
    int lmax, mmax;
    std::unordered_map<int,                 // l
        std::unordered_map<int,             // m
            std::vector<Population>>>       // populations
                populations;
    std::unordered_map<int, std::vector<std::pair<double, double>>> phaseshifts;

    bool Validate(const nlohmann::json& input);
    void Load(const nlohmann::json& input);
    void Execute();
    void CalculatePhaseshifts();

    void AngularlyResolved();
    void EnergyResolved();
    void ReadPopulations();
};




#endif
