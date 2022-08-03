#include "objects/observables/pulse_obs.h"
#include "common/tdse/simulation.h"
#include "common/file_io/io_factory.h"
#include "common/utility/logger.h"
#include <sstream>
#include <iomanip>

using namespace tdse;
using Maths = maths::Factory;
using namespace io;

PulseObservable::PulseObservable() {}         

void PulseObservable::Startup(int it) {
    // dump the whole field into a text file
    if (_outputFilename.length() > 0) {
        _txtFile = ASCII(io::Factory::OpenASCII(_outputFilename, 'w'));
    }

    // set up stringstream for formating and grab some shortcut to values we need
    std::stringstream ss;
    ss << std::setprecision(8) << std::scientific;
    double dt = Simulation::GetTimestep();
    int NT = Simulation::GetNumTimeSteps();
    auto& pulses = Simulation::GetPulses();
    // auto& A_field_x = Simulation::GetField(DimIndex::X);
    // auto& A_field_y = Simulation::GetField(DimIndex::Y);
    // auto& A_field_z = Simulation::GetField(DimIndex::Z);
    
    std::vector<Vec3> A_field(NT, {0,0,0});
    std::vector<Vec3> E_field(NT, {0,0,0});

    for (auto& pulse : pulses) {
        for (int i = 0; i < NT; i++) {
            A_field[i] += pulse->A(dt*i);
            E_field[i] += pulse->E(dt*i);
        }
    }

    // loop though the entire pulse
    for (int i = 0; i < NT; i++) {

        ss.str("");
        ss << dt*i;
        // A-field
        /*
        if (Simulation::Polarization()[DimIndex::X])
            ss << "\t" << A_field[i].x;
        else 
            ss << "\t" << 0.0;
        if (Simulation::Polarization()[DimIndex::Y])
            ss << "\t" << A_field[i].y;
        else 
            ss << "\t" << 0.0;
        if (Simulation::Polarization()[DimIndex::Z])
            ss << "\t" << A_field[i].z;
        else 
            ss << "\t" << 0.0;
        */

        ss << "\t" << A_field[i].x;
        ss << "\t" << A_field[i].y;
        ss << "\t" << A_field[i].z;
        ss << "\t" << E_field[i].x;
        ss << "\t" << E_field[i].y;
        ss << "\t" << E_field[i].z;

        ss << "\n";
        _txtFile->Write(ss.str());
    }

    // done
    _txtFile = nullptr;
}
void PulseObservable::Shutdown() {}
void PulseObservable::Compute(int it) {}
void PulseObservable::Flush() {}


std::string PulseObservable::GetName() {
    return "pulse";
}
tdse::Observable::Ptr_t PulseObservable::Create(const nlohmann::json& observable) {
    auto pl_obs = new PulseObservable();
    pl_obs->_outputFilename = observable["filename"].get<std::string>();
    return tdse::Observable::Ptr_t(pl_obs);
}
bool PulseObservable::Validate(const nlohmann::json& observable) {
    // an output file is required
    if (!observable.contains("filename")) {
        LOG_CRITICAL("\"pulse\" observable must contain string entry: filename");
        return false;
    }
    return true;
}

