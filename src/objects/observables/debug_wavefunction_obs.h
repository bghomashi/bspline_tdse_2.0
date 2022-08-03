#ifndef __DEBUG_WAVEFUNCTION_H__
#define __DEBUG_WAVEFUNCTION_H__

#include "common/objects/observable.h"
#include "common/file_io/io_ascii.h"
#include <string>

class DebugWavefunctionObservable : public tdse::Observable::Register<DebugWavefunctionObservable> {
    io::ASCII _txtFile;
    std::string _outputFilename;

    int _numGrid;
    std::vector<double> _grid;
public:
    DebugWavefunctionObservable();

    void SetNumGrid(int numGrid);

    void Flush();
    void Startup(int it);
    void Shutdown();
    void Compute(int it);

    static std::string GetName();
    static Observable::Ptr_t Create(const nlohmann::json& observable);
    static bool Validate(const nlohmann::json& observable);
};

#endif