#include "objects/observables/debug_wavefunction_obs.h"
#include "common/utility/logger.h"
#include "common/tdse/simulation.h"
#include "common/file_io/io_factory.h"
#include "common/maths/math_common.h"
#include "common/bspline/bspline.h"
#include <iomanip>

using namespace tdse;
using namespace maths;

DebugWavefunctionObservable::DebugWavefunctionObservable() : _numGrid(0) {}

void DebugWavefunctionObservable::SetNumGrid(int numGrid) {
    _numGrid = numGrid;
}

void DebugWavefunctionObservable::Startup(int start_it) {
    _grid.resize(_numGrid);

    double dx = (bspline::Basis::GetXmax() - bspline::Basis::GetXmin())/(_numGrid - 1);
    for (int i = 0; i < _numGrid; i++)
        _grid[i] = bspline::Basis::GetXmin() + i*dx;

}
void DebugWavefunctionObservable::Shutdown() {
    _txtFile = io::Factory::OpenASCII("wavefunction_final.txt", 'w');

    std::vector<complex> psi(Simulation::GetDOF()), n_block_coeff;
    std::vector<complex> out(_grid.size(), 0), temp(_grid.size());
    std::stringstream ss;
    
    Simulation::GetPsi()->CopyTo(psi);

    int N = bspline::Basis::GetNumBSplines();

    for (int i = 0; i < Simulation::GetDOF(); i += N) {
        n_block_coeff = std::vector<complex>(psi.begin() + i, psi.begin() + (i+N));             // copy n_block coeffs
        temp = bspline::Basis::FunctionEvaluate(_grid, n_block_coeff);                                        // evaluate on grid
        
        for (int i = 0; i < _grid.size(); i++)                                                      // sum the square at all the grid points                  
            out[i] += std::abs(temp[i])*std::abs(temp[i]);
    }
    ss << std::setprecision(8) << std::scientific;
    for (int i = 0; i < _numGrid; i++) {            // n-quantum number
        ss.str("");
        ss << _grid[i] << "\t" 
           << std::real(out[i]) << "\n";
        _txtFile->Write(ss.str());
    }


    _txtFile = nullptr;
}
void DebugWavefunctionObservable::Compute(int it) {
    _txtFile = io::Factory::OpenASCII("wavefunction_"+std::to_string(it)+".txt", 'w');

    std::vector<complex> psi(Simulation::GetDOF()), n_block_coeff;
    std::vector<complex> out(_grid.size(), 0), temp(_grid.size());
    std::stringstream ss;
    
    Simulation::GetPsi()->CopyTo(psi);

    int N = bspline::Basis::GetNumBSplines();

    for (int i = 0; i < Simulation::GetDOF(); i += N) {
        n_block_coeff = std::vector<complex>(psi.begin() + i, psi.begin() + (i+N));             // copy n_block coeffs
        temp = bspline::Basis::FunctionEvaluate(_grid, n_block_coeff);                                        // evaluate on grid
        
        for (int i = 0; i < _grid.size(); i++)                                                      // sum the square at all the grid points                  
            out[i] += std::abs(temp[i])*std::abs(temp[i]);
    }

    for (int i = 0; i < _numGrid; i++) {            // n-quantum number
        ss.str("");
        ss << _grid[i] << "\t" 
           << std::real(out[i]) << "\n";
        _txtFile->Write(ss.str());
    }


    _txtFile = nullptr;
}


void DebugWavefunctionObservable::Flush() {
}

std::string DebugWavefunctionObservable::GetName() {
    return "debug_wavefunction";
}
Observable::Ptr_t DebugWavefunctionObservable::Create(const nlohmann::json& observable) {
    auto wf_obs = new DebugWavefunctionObservable();
    if (observable.contains("compute_period")) 
        wf_obs->_computePeriod = observable["compute_period"].get<int>();
    wf_obs->_numGrid = observable["grid_points"].get<int>();
    return Observable::Ptr_t(wf_obs);
}
bool DebugWavefunctionObservable::Validate(const nlohmann::json& observable) {
    if (observable.contains("compute_period") && !observable["compute_period"].is_number()) {
        LOG_CRITICAL("Optional entry \"compute_period\" must be a number.");
        return false;
    }
    if (!(observable.contains("grid_points") && observable["grid_points"].is_number())) {
        LOG_CRITICAL("\"debug_wavefunction\" observable must contain number: grid_points");
        return false;
    }
    return true;
}