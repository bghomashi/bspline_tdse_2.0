#pragma once

#include "common/objects/observable.h"
#include "common/maths/math_common.h"
#include "common/file_io/io_ascii.h"

class IonizationObservable : public tdse::Observable::Register<IonizationObservable> {
    std::string _outputFilename;
    
    void IonizationCentral();
    void IonizationAxial();
    void BuildOverlap(maths::Matrix S);
public:
    IonizationObservable();

    void Startup(int it);
    void Shutdown();
    void Compute(int it);
    void Flush();
    

    static std::string GetName();
    static Observable::Ptr_t Create(const nlohmann::json& observable);
    static bool Validate(const nlohmann::json& observable);
};
