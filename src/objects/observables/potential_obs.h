#ifndef __POTENTIAL_OBS_H__
#define __POTENTIAL_OBS_H__

#include "common/objects/observable.h"
#include "common/file_io/io_ascii.h"

class PotentialObservable : public tdse::Observable::Register<PotentialObservable> {
protected:
    io::ASCII _txtFile;
    int _numGrid;

    std::string _outputFilename;
public:
    PotentialObservable();         
    
    void Flush();
    void Startup(int it);
    void Shutdown();
    void Compute(int it);

    static std::string GetName();
    static Observable::Ptr_t Create(const nlohmann::json& observable);
    static bool Validate(const nlohmann::json& observable);
};

#endif