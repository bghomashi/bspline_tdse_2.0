#ifndef __BASIS_OBS_H__
#define __BASIS_OBS_H__

#include "common/objects/observable.h"
#include "common/file_io/io_ascii.h"
#include <vector>

class BasisObservable : public tdse::Observable::Register<BasisObservable> {
protected:
    std::string _outputFilename;
    io::ASCII _txtFile;
    int _numGrid, _to, _from;
    std::vector<double> _grid;
public:
    BasisObservable();
        
    void Flush();
    void Startup(int it);
    void Shutdown();
    void Compute(int it);

    static std::string GetName();
    static Observable::Ptr_t Create(const nlohmann::json& observable);
    static bool Validate(const nlohmann::json& observable);
};

#endif