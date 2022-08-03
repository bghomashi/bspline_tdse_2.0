#ifndef __DEBUG_EIGENSTATES_H__
#define __DEBUG_EIGENSTATES_H__

#include "common/file_io/io_ascii.h"
#include "common/objects/observable.h"
#include <string>

class DebugEigenstatesObservable : public tdse::Observable::Register<DebugEigenstatesObservable> {
    io::ASCII _txt_file;

    int _numGrid, _nmax;
    std::vector<double> _grid;
public:
    DebugEigenstatesObservable();

    void SetNumGrid(int numGrid);
    void SetNMax(int nmax);
    void Startup(int it);
    void Shutdown();
    void Compute(int it);

    static std::string GetName();
    static Observable::Ptr_t Create(const nlohmann::json& observable);
    static bool Validate(const nlohmann::json& observable);
};

#endif