#include "objects/observables/wavefunction_obs.h"
#include "common/maths/math_factory.h"
#include "common/file_io/io_factory.h"
#include "common/tdse/simulation.h"
#include "common/utility/logger.h"
#include "common/bspline/bspline.h"
#include "common/system_state/system_state.h"

using namespace tdse;
using Math = maths::Factory;
using namespace maths;

WavefunctionObservable::WavefunctionObservable() : _numGrid(0) {}

void WavefunctionObservable::Startup(int start_it) {
    _psi = Simulation::GetPsi();
    _psiGrid = Math::CreateVector(_numGrid);
    _grid.resize(_numGrid);

    double dx = (bspline::Basis::GetXmax() - bspline::Basis::GetXmin())/(_numGrid - 1);
    for (int i = 0; i < _numGrid; i++)
        _grid[i] = bspline::Basis::GetXmin() + i*dx;

    _hdf5 = io::Factory::OpenHDF5(_outputFilename, 'w');
}
void WavefunctionObservable::Shutdown() {
    _psi = nullptr;
    _psiGrid = nullptr;
    _hdf5 = nullptr;
}
void WavefunctionObservable::Compute(int it) {
    std::stringstream ss;
    double t = Simulation::GetTime();
    _psi->Transform(_psiGrid, [=](const std::vector<complex>& coeff) {
        auto v = bspline::Basis::FunctionEvaluate(_grid, coeff);
        return v;
    });
    
    _hdf5->PushGroup("vectors");
    for (int l = 0; l < SystemState::GetBasisLmax(); l++) {            // l-quantum number
        for (int i = 0; i < _numGrid; i++) {            // n-quantum number
            ss.str("");                                 // clear string stream
            ss << "psi[" << it << "]";

            _hdf5->WriteVector(ss.str().c_str(), _psiGrid);
            _hdf5->WriteAttribute(_psiGrid, "time", t);
        }
    }    
    _hdf5->PopGroup();
}


std::string WavefunctionObservable::GetName() {
    return "wavefunction";
}
Observable::Ptr_t WavefunctionObservable::Create(const nlohmann::json& observable) {
    auto wf_obs = new WavefunctionObservable();
    if (observable.contains("compute_period")) 
        wf_obs->_computePeriod = observable["compute_period"].get<int>();
    wf_obs->_outputFilename = observable["filename"].get<std::string>();
    wf_obs->_numGrid = observable["grid_points"].get<int>();

    return Observable::Ptr_t(wf_obs);
}
bool WavefunctionObservable::Validate(const nlohmann::json& observable) {
    if (observable.contains("compute_period") && !observable["compute_period"].is_number()) {
        LOG_CRITICAL("Optional entry \"compute_period\" must be a number.");
        return false;
    }
    if (!(observable.contains("grid_points") && observable["grid_points"].is_number())) {
        LOG_CRITICAL("\"wavefunction\" observable must contain number: grid_points");
        return false;
    }
    if (!observable.contains("filename")) {
        LOG_CRITICAL("\"wavefunction\" observable must contain string entry: filename");
        return false;
    }
    return true;
}
