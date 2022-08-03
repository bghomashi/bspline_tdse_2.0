#include "objects/observables/debug_eigenstates.h"
#include "common/utility/logger.h"
#include "common/maths/math_common.h"
#include "common/maths/math_factory.h"
#include "common/file_io/io_factory.h"
#include "common/bspline/bspline.h"
#include "common/tdse/simulation.h"
#include <iomanip>

using namespace tdse;
using Math = maths::Factory;
using namespace maths;

DebugEigenstatesObservable::DebugEigenstatesObservable() : _numGrid(0) {}

void DebugEigenstatesObservable::SetNumGrid(int numGrid) {
    _numGrid = numGrid;
}

void DebugEigenstatesObservable::Startup(int start_it) {
    _grid.resize(_numGrid);

    double dx = (bspline::Basis::GetXmax() - bspline::Basis::GetXmin())/(_numGrid - 1);
    for (int i = 0; i < _numGrid; i++)
        _grid[i] = bspline::Basis::GetXmin() + i*dx;

    std::stringstream name_ss, ss;
    std::vector<complex> temp(_grid.size()), out(_grid.size());
    int N = bspline::Basis::GetNumBSplines();
    Vector eigen = Math::CreateVector(N);        // vector from hdf5 file

    auto hdf5 = io::Factory::OpenHDF5(Simulation::GetInitialStateFile(), 'r');
    hdf5->PushGroup("vectors");

    ss << std::setprecision(8) << std::scientific;
    for (int n = 1; n <= _nmax; n++) {
        for (int l = 0; l < n; l++) {
            name_ss.str("");                                         // clear string stream
            name_ss << "(" << n << ", " << l << ")";     // name of state

            _txt_file = io::Factory::OpenASCII("eigenstate_"+name_ss.str()+".txt", 'w');

            hdf5->ReadVector(name_ss.str().c_str(), eigen);           // read in the state
            eigen->CopyTo(temp);
            out = bspline::Basis::FunctionEvaluate(_grid, temp);  

            for (int i = 0; i < _numGrid; i++) {            // n-quantum number
                ss.str("");
                ss << _grid[i] << "\t" <<
                      std::real(out[i]) << "\n";
                _txt_file->Write(ss.str());
            }                                      // evaluate on grid
        }
    }
    hdf5->PopGroup();

    hdf5 = nullptr;
    eigen = nullptr;
    _txt_file = nullptr;
}

void DebugEigenstatesObservable::SetNMax(int nmax) {
    _nmax = nmax;
}
void DebugEigenstatesObservable::Compute(int it) {}
void DebugEigenstatesObservable::Shutdown() {}



std::string DebugEigenstatesObservable::GetName() {
    return "debug_eigenstates";
}
Observable::Ptr_t DebugEigenstatesObservable::Create(const nlohmann::json& observable) {
    auto eigen_obs = new DebugEigenstatesObservable();
    eigen_obs->_numGrid = observable["grid_points"].get<int>();
    eigen_obs->_nmax = observable["nmax"].get<int>();
    return Observable::Ptr_t(eigen_obs);
}
bool DebugEigenstatesObservable::Validate(const nlohmann::json& observable) {
    if (!(observable.contains("nmax") && observable["nmax"].is_number())) {
        LOG_CRITICAL("\"debug_eigenstates\" observable must contain number entry: nmax");
        return false;
    }
    if (!(observable.contains("grid_points") && observable["grid_points"].is_number())) {
        LOG_CRITICAL("\"debug_eigenstates\" observable must contain number entry: grid_points");
        return false;
    }
    return true;
}