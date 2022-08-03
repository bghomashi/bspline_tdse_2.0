#ifndef __PULSE_OBS_H__
#define __PULSE_OBS_H__

#include "common/objects/observable.h"
#include "common/file_io/io_ascii.h"

class PulseObservable : public tdse::Observable::Register<PulseObservable> {
protected:
    io::ASCII _txtFile;
    std::string _outputFilename;
public:
    PulseObservable();         
    
    void Flush();
    void Startup(int it);
    void Shutdown();
    void Compute(int it);

    static std::string GetName();
    static tdse::Observable::Ptr_t Create(const nlohmann::json& observable);
    static bool Validate(const nlohmann::json& observable);
};

#endif