#include "objects/observables/density_obs.h"
#include "common/utility/logger.h"
#include "common/utility/index_manip.h"
#include "common/utility/spherical_harmonics.h"
#include "common/tdse/simulation.h"
#include "common/maths/math_factory.h"
#include "common/file_io/io_factory.h"
#include "common/bspline/bspline.h"
#include "common/system_state/system_state.h"
#include <iomanip>
#include <cmath>

using namespace tdse;
using Maths = maths::Factory;
using namespace maths;

DensityObservable::DensityObservable() : _numGrid(0) {}

void DensityObservable::SetNumGrid(int numGrid) {
    _numGrid = numGrid;
}

void DensityObservable::Startup(int start_it) {
    // same grid in all dimensions
    _grid.resize(_numGrid);

    // make the grid symmetric [-max, max]
    double dx = 2.*bspline::Basis::GetXmax()/(_numGrid - 1);
    for (int i = 0; i < _numGrid; i++)
        _grid[i] = -bspline::Basis::GetXmax() + i*dx;

}
void DensityObservable::Shutdown() {
    _txt_file = nullptr;
}
void DensityObservable::Compute(int it) {
    _txt_file = io::Factory::OpenASCII("density_"+std::to_string(it)+".txt", 'w');

    std::vector<maths::complex> psi(Simulation::GetDOF()), n_block_coeff;
    std::vector<maths::complex> temp(_grid.size());
    std::stringstream ss;
    
    int N = bspline::Basis::GetNumBSplines();
    int lmax = SystemState::GetBasisLmax();
    auto& Ms = Simulation::GetMs();
    auto& MRows = Simulation::GetMRows();
    

    Simulation::GetPsi()->CopyTo(psi);                                                                   // copy entire wf into vector

    for (int i = 0; i < _numGrid; i++) {
        for (int j = 0; j < _numGrid; j++) {
            for (int k = 0; k < _numGrid; k++) {
                double x = _grid[i];
                double y = _grid[j];
                double z = _grid[k];

                double r = std::sqrt(x*x + y*y + z*z);
                double theta = (x == 0 && y == 0 && z == 0 ? 0.0 : std::atan2(std::sqrt(x*x + y*y), z));
                double phi = (x == 0 && y == 0 ? 0.0 : std::atan2(y, x));

                // if r == 0 : r = epsilon ? to avoid divide by zero?

                complex amplitude = 0.0;
                for (int m : Ms) {                                                                              // for each m
                    for (int l = std::abs(m); l <= lmax; l++) {                                                 // for each l
                        int start = RowFrom(0, m, l, N, Ms, MRows);                                             // first index of this block
                        n_block_coeff = std::vector<complex>(psi.begin() + start, psi.begin() + (start+N)); // copy n_block coeffs
            
                        amplitude += bspline::Basis::FunctionEvaluate(r, n_block_coeff)*Ylm(l, m, theta, phi) / r; // evaluate on grid this chuck on the grid
                    }
                }

                ss.str("");
                ss << x << "\t";
                ss << y << "\t";
                ss << z << "\t";
                ss << std::abs(amplitude)*std::abs(amplitude) << "\t";
                ss << std::real(amplitude) << "\t";
                ss << std::imag(amplitude) << "\n";
                _txt_file->Write(ss.str()); 

            }            
        }  
    }       

    _txt_file = nullptr;
}


std::string DensityObservable::GetName() {
    return "density";
}
Observable::Ptr_t DensityObservable::Create(const nlohmann::json& observable) {
    auto density_obs = new DensityObservable();
    if (observable.contains("compute_period")) 
        density_obs->_computePeriod = observable["compute_period"];
    density_obs->_numGrid = observable["grid_points"];
    return Observable::Ptr_t(density_obs);
}
bool DensityObservable::Validate(const nlohmann::json& observable) {
    if (observable.contains("compute_period") && !observable["compute_period"].is_number()) {
        LOG_CRITICAL("Optional entry \"compute_period\" must be a number.");
        return false;
    }
     if (!(observable.contains("grid_points") && observable["grid_points"].is_number())) {
        LOG_CRITICAL("\"density\" observable must contain number entry: grid_points");
        return false;
    }
    return true;
}