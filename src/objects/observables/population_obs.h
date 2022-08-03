#pragma once

#include "common/objects/observable.h"
#include "common/maths/math_common.h"
#include "common/file_io/io_ascii.h"

class PopulationObservable : public tdse::Observable::Register<PopulationObservable> {
    std::string _outputFilename;
    
    void PopulationsCentral();
    void PopulationsAxial();
    void BuildOverlap(maths::Matrix S);
public:
    PopulationObservable();

    void Startup(int it);
    void Shutdown();
    void Compute(int it);
    void Flush();
    

    static std::string GetName();
    static Observable::Ptr_t Create(const nlohmann::json& observable);
    static bool Validate(const nlohmann::json& observable);
};
