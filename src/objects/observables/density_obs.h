#ifndef __DENSITY_OBS_H__
#define __DENSITY_OBS_H__

#include <string>
#include "common/objects/observable.h"
#include "common/file_io/io_ascii.h"

// this while output the 3D density in cartesian coordinates. Use responsibly
class DensityObservable : public tdse::Observable::Register<DensityObservable> {
    io::ASCII _txt_file;

    int _numGrid;
    std::vector<double> _grid;
public:
    DensityObservable();

    void SetNumGrid(int numGrid);

    void Startup(int it);
    void Shutdown();
    void Compute(int it);

    static std::string GetName();
    static Observable::Ptr_t Create(const nlohmann::json& observable);
    static bool Validate(const nlohmann::json& observable);
};

#endif