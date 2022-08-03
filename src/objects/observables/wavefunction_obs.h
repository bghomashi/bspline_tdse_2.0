#ifndef __WAVEFUNCTION_OBS_H__
#define __WAVEFUNCTION_OBS_H__

#include "common/objects/observable.h"
#include "common/file_io/io_hdf5.h"
#include "common/maths/math_common.h"
#include <string>

class WavefunctionObservable : public tdse::Observable::Register<WavefunctionObservable> {
    maths::Vector _psi;            // shortcut to wavefunction
    maths::Vector _psiGrid;        // holds wavefunction in real-space 
    io::HDF5 _hdf5;
    std::string _outputFilename;

    int _numGrid;
    std::vector<double> _grid;
public:
    WavefunctionObservable();

    void Startup(int it);
    void Shutdown();
    void Compute(int it);

    static std::string GetName();
    static tdse::Observable::Ptr_t Create(const nlohmann::json& observable);
    static bool Validate(const nlohmann::json& observable);
    
};

#endif