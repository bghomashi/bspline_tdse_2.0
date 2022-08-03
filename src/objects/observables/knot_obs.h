#ifndef __KNOTS_OBS_H__
#define __KNOTS_OBS_H__

#include "common/objects/observable.h"
#include "common/file_io/io_ascii.h"

class KnotObservable : public tdse::Observable::Register<KnotObservable> {
protected:
    std::string _outputFilename;
    io::ASCII _txtFile;
public:
    KnotObservable();
        
    void Flush();
    void Startup(int it);
    void Shutdown();
    void Compute(int it);

    static std::string GetName();
    static Observable::Ptr_t Create(const nlohmann::json& observable);
    static bool Validate(const nlohmann::json& observable);
};

#endif