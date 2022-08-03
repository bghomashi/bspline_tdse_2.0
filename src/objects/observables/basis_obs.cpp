#include "objects/observables/basis_obs.h"
#include "common/maths/math_common.h"
#include "common/file_io/io_factory.h"
#include "common/utility/logger.h"
#include "common/bspline/bspline.h"
#include <sstream>
#include <iomanip>
#include <limits.h>

using namespace tdse;
using namespace maths;

BasisObservable::BasisObservable() : _from(0), _to(INT_MAX) {}         

void BasisObservable::Startup(int it) {
    // dump the whole field into a text file
    if (_outputFilename.length() > 0)
        _txtFile = io::Factory::OpenASCII(_outputFilename, 'w');

    int num_bsplines = bspline::Basis::GetNumBSplines();
    _to = std::min(num_bsplines-1, _to);
    std::vector<double> grid(_numGrid);
    std::vector<std::vector<complex>> splines(_to - _from + 1);

    double dx = (bspline::Basis::GetXmax() - bspline::Basis::GetXmin())/(_numGrid - 1);
    for (int i = 0; i < _numGrid; i++)
        grid[i] = bspline::Basis::GetXmin() + i*dx;
        
    for (int bs = 0; bs < splines.size(); bs++)
        splines[bs] = bspline::Basis::GetBSpline(grid, bs + _from);

    // set up stringstream for formating and grab some shortcut to values we need
    std::stringstream ss;
    ss << std::setprecision(8) << std::scientific;
    
    // loop though the entire pulse
    for (int i = 0; i < _numGrid; i++) {
        ss.str("");
        ss << grid[i];
        for (auto& bs : splines) {
            ss << "\t" << std::real(bs[i]);
            ss << "\t" << std::imag(bs[i]);
        }
        ss << "\n";
        _txtFile->Write(ss.str());
    }

    // done
    _txtFile = nullptr;
}
void BasisObservable::Shutdown() {}
void BasisObservable::Compute(int it) {}
void BasisObservable::Flush() {}


std::string BasisObservable::GetName() {
    return "basis";
}
Observable::Ptr_t BasisObservable::Create(const nlohmann::json& observable) {
    auto basis_obs = new BasisObservable();
    basis_obs->_numGrid = observable["grid_points"].get<int>();
    basis_obs->_outputFilename = observable["filename"].get<std::string>();
    if (observable.contains("from"))
        basis_obs->_from = observable["from"].get<int>();
    if (observable.contains("to"))
        basis_obs->_to = observable["to"].get<int>();
    return Observable::Ptr_t(basis_obs);
}
bool BasisObservable::Validate(const nlohmann::json& observable) {
    if (!(observable.contains("grid_points") && observable["grid_points"].is_number())) {
        LOG_CRITICAL("\"basis\" observable must contain number entry: grid_points");
        return false;
    }
    if (observable.contains("from") && !observable["from"].is_number()) {
        LOG_CRITICAL("\"basis\" observable must contain number entry: from");
        return false;
    }
    if (observable.contains("to") && !observable["to"].is_number()) {
        LOG_CRITICAL("\"basis\" observable must contain number entry: to");
        return false;
    }
    // an output file is required
    if (!observable.contains("filename")) {
        LOG_CRITICAL("\"basis\" observable must contain string entry: filename");
        return false;
    }
    return true;
}