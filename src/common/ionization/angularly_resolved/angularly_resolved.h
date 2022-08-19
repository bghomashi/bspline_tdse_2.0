#ifndef __ANGULARLY_RESOLVED_H__
#define __ANGULARLY_RESOLVED_H__

#include <string>
#include "common/ionization/ionization_task/ionization_task.h"
#include "common/utility/json.hpp"

class AngularlyResolved : public IonizationTask::Register<AngularlyResolved> {
    std::string _output_filename;
    double _emin, _emax, _estep;
    double _xmin, _xmax, _xstep;
    double _zmin, _zmax, _zstep;
public:
    static std::string GetName();
    static bool Validate(const nlohmann::json& angularly_resolved);
    static IonizationTask::Ptr_t Create(const nlohmann::json& angularly_resolved);

    void Execute();
};


#endif