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
    double dx = (_xmax - _xmin)/(_numGrid - 1);
    for (int i = 0; i < _numGrid; i++)
        _grid[i] = _xmin + i*dx;

}
void DensityObservable::Shutdown() {
    _bin_file = nullptr;
}
void DensityObservable::Compute(int it) {
    _bin_file = io::Factory::OpenBinary("density_"+std::to_string(it)+".bin", 'w');

    std::vector<maths::complex> psi, n_block_coeff;
    std::vector<maths::complex> temp(_grid.size());
    std::stringstream ss;
    
    int N = bspline::Basis::GetNumBSplines();
    int lmax = SystemState::GetBasisLmax();
    auto& Ms = Simulation::GetMs();
    auto& MRows = Simulation::GetMRows();
    double dx = (_xmax - _xmin)/(_numGrid - 1);
    
    Simulation::GetPsi()->CopyTo(psi);                                                                   // copy entire wf into vector

    if (psi.size() > 0) {
        _bin_file->Write(&_numGrid, sizeof(int)); 
        _bin_file->Write(&dx, sizeof(double)); 
        _bin_file->Write(&_grid[0], sizeof(double)); 
        _bin_file->Write(&_grid[_grid.size()-1], sizeof(double)); 

        for (int i = 0; i < _numGrid; i++) {
            for (int j = 0; j < _numGrid; j++) {
                // for (int k = 0; k < _numGrid; k++)
                {
                    double x = _grid[i];
                    double y = _grid[j];
                    // double z = _grid[k];
                    double z = 0;

                    double r = std::sqrt(x*x + y*y + z*z);
                    double theta = (x == 0 && y == 0 && z == 0 ? 0.0 : std::atan2(std::sqrt(x*x + y*y), z));
                    double phi = (x == 0 && y == 0 ? 0.0 : std::atan2(y, x));

                    if (r == 0) 
                        r = FLT_MIN;

                    complex amplitude = 0.0;
                    for (int m : Ms) {                                                                              // for each m
                        for (int l = std::abs(m); l <= lmax; l++) {                                                 // for each l
                            int start = RowFrom(0, l, m, N, Ms, MRows);                                             // first index of this block
                            n_block_coeff = std::vector<complex>(psi.begin() + start, psi.begin() + (start+N)); // copy n_block coeffs
                
                            amplitude += bspline::Basis::FunctionEvaluate(r, n_block_coeff)*Ylm(l, m, theta, phi) / r; // evaluate on grid this chuck on the grid
                            if (std::isnan(std::real(amplitude)) || std::isnan(std::imag(amplitude))) {
                                std::cout << "isnan: " << i << " " << j << std::endl;
                                std::cout << "l " << l << " m " << m << std::endl;
                                std::cout << "bspline " << bspline::Basis::FunctionEvaluate(r, n_block_coeff) << std::endl;
                                std::cout << "RowFrom(0, m, l, N, Ms, MRows); " << RowFrom(0, l, m, N, Ms, MRows) << std::endl;
                                exit(0);
                            }
                            // if (std::abs(amplitude) > 0) {
                            //     std::cout << i << " " << j << std::endl;
                            // }
                        }
                    }

                    _bin_file->Write(&amplitude, sizeof(complex)); 

                }            
            }  
        }
        _bin_file->Flush();
        _bin_file = nullptr;
    }
}


std::string DensityObservable::GetName() {
    return "density";
}
Observable::Ptr_t DensityObservable::Create(const nlohmann::json& observable) {
    auto density_obs = new DensityObservable();
    if (observable.contains("compute_period")) 
        density_obs->_computePeriod = observable["compute_period"];
    density_obs->_numGrid = observable["grid_points"];
    if (observable.contains("xmin")) 
        density_obs->_xmin = observable["xmin"];
    else 
        density_obs->_xmin = -bspline::Basis::GetXmax();
    if (observable.contains("xmax")) 
        density_obs->_xmax = observable["xmax"];
    else 
        density_obs->_xmin = bspline::Basis::GetXmax();
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