#include "objects/observables/knot_obs.h"
#include "common/tdse/simulation.h"
#include "common/file_io/io_factory.h"
#include "common/utility/logger.h"
#include "common/bspline/bspline.h"
#include <sstream>
#include <iomanip>
#include <limits.h>

using namespace tdse;

KnotObservable::KnotObservable() {}         

void KnotObservable::Startup(int it) {
    // dump the whole field into a text file
    if (_outputFilename.length() > 0) {
        _txtFile = io::ASCII(io::Factory::OpenASCII(_outputFilename, 'w'));
    }

    auto& grid = bspline::Basis::GetGrid();
    
    // set up stringstream for formating and grab some shortcut to values we need
    std::stringstream ss;
    ss << std::setprecision(8) << std::scientific;
    
    // loop though the entire pulse
    for (int i = 0; i < grid.size(); i++) {
        ss.str("");
        ss << i << "\t"
           << grid[i] << "\n";
        _txtFile->Write(ss.str());
    }

    // done
    _txtFile = nullptr;
}



void KnotObservable::Flush() {}
void KnotObservable::Shutdown() {}
void KnotObservable::Compute(int it) {}


std::string KnotObservable::GetName() {
    return "knots";
}
Observable::Ptr_t KnotObservable::Create(const nlohmann::json& observable) {
    auto knots_obs = new KnotObservable();

    knots_obs->_outputFilename = observable["filename"];
    
    return Observable::Ptr_t(knots_obs);
}
bool KnotObservable::Validate(const nlohmann::json& observable) {
    // an output file is required
    if (!observable.contains("filename")) {
        LOG_CRITICAL("\"knots\" observable must contain string entry: filename");
        return false;
    }
    return true;
}

