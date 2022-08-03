#include "objects/observables/potential_obs.h"
#include "common/file_io/io_factory.h"
#include "common/tdse/simulation.h"
#include "common/utility/logger.h"
#include "common/bspline/bspline.h"
#include <sstream>
#include <iomanip>

using namespace tdse;

PotentialObservable::PotentialObservable() {}         

void PotentialObservable::Startup(int it) {
    // dump the whole field into a text file
    if (_outputFilename.length() > 0) {
        _txtFile = io::ASCII(io::Factory::OpenASCII(_outputFilename, 'w'));
    }

    // set up stringstream for formating and grab some shortcut to values we need
    std::stringstream ss;
    ss << std::setprecision(8) << std::scientific;
    auto& potentials = Simulation::GetPotentials();
    double xmin = bspline::Basis::GetXmin();
    double xmax = bspline::Basis::GetXmax();
    double dx = (xmax - xmin)/(_numGrid - 1.);
    
    // this assumes all potentials are central!
    // loop though the entire pulse
    for (int i = 0; i < _numGrid; i++) {
        double r = xmin + dx*i;
        double p = 0;
        for (auto& pot : potentials)
            p += (*pot)(r,0,0);

        ss.str("");
        ss << r << "\t" << p << "\n";
        _txtFile->Write(ss.str());
    }

    // done
    _txtFile = nullptr;
}

void PotentialObservable::Shutdown() {}
void PotentialObservable::Compute(int it) {}
void PotentialObservable::Flush() {}



std::string PotentialObservable::GetName() {
    return "potential";
}
Observable::Ptr_t PotentialObservable::Create(const nlohmann::json& observable) {
    auto pot_obs = new PotentialObservable();
    pot_obs->_numGrid = observable["grid_points"].get<int>();
    pot_obs->_outputFilename = observable["filename"].get<std::string>();
    return Observable::Ptr_t(pot_obs);
}
bool PotentialObservable::Validate(const nlohmann::json& observable) {
    if (!(observable.contains("grid_points") && observable["grid_points"].is_number())) {
        LOG_CRITICAL("\"potential\" observable must contain number entry: grid_points");
        return false;
    }
    // an output file is required
    if (!observable.contains("filename")) {
        LOG_CRITICAL("\"potential\" observable must contain string entry: filename");
        return false;
    }
    return true;
}