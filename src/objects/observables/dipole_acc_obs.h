#ifndef __DIPOLE_ACC_OBS_H__
#define __DIPOLE_ACC_OBS_H__

#include "common/objects/observable.h"
#include "common/maths/math_common.h"
#include "common/file_io/io_ascii.h"

class DipoleAccObservable : public tdse::Observable::Register<DipoleAccObservable> {
    maths::Matrix _gradPot[maths::DimIndex::NUM];
    maths::Vector _psi;            // shortcut to wavefunction
    maths::Vector _psiTemp;       // just from MatMult output

    std::string _outputFilename;
    io::ASCII _txtFile;
public:
    DipoleAccObservable();

    int MemoryAlloced() const;
    void Flush();
    void Startup(int it);
    void Shutdown();
    void Compute(int it);

    static std::string GetName();
    static Observable::Ptr_t Create(const nlohmann::json& observable);
    static bool Validate(const nlohmann::json& observable);
};

#endif