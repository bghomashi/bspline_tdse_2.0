#pragma once

#include "common/objects/observable.h"
#include "common/maths/math_common.h"
#include "common/file_io/io_ascii.h"

class NormObservable : public tdse::Observable::Register<NormObservable> {
    maths::Matrix _S;
    maths::Vector _psi;            // shortcut to wavefunction
    maths::Vector _psi_temp;

    std::string _outputFilename;
    io::ASCII _txtFile;
public:
    NormObservable();

    void Startup(int it);
    void Shutdown();
    void Compute(int it);
    void Flush();

    static std::string GetName();
    static Observable::Ptr_t Create(const nlohmann::json& observable);
    static bool Validate(const nlohmann::json& observable);
};
