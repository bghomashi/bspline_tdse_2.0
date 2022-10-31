#include "objects/observables/norm_obs.h"
#include "common/utility/logger.h"
#include "common/utility/file_exists.h"
#include "common/utility/banded_matrix.h"
#include "common/tdse/simulation.h"
#include "common/maths/math_factory.h"
#include "common/maths/math_algebra.h"
#include "common/file_io/io_factory.h"
#include "common/bspline/bspline.h"
#include <iostream>
#include <sstream>

using namespace tdse;
using Math = maths::Factory;
using namespace maths;

NormObservable::NormObservable() {}

void NormObservable::Startup(int start_it) {
    int order = bspline::Basis::GetOrder();
    int N = bspline::Basis::GetNumBSplines();
    int dof = Simulation::GetDOF();

    _S = Math::CreateMatrix(dof, dof, 2*order-1);
    _psi = Simulation::GetPsi();
    _psi_temp = Math::CreateVector(dof);

    // --------------- Fill overlap matrix
    BandedMatrix overlapStore(N, 2*order-1);
    for (int r = 0; r < N; r++) {
        for (int c = r; c < std::min(N, r+order); c++) {
            overlapStore(r,c) = overlapStore(c,r) = bspline::Basis::Integrate(r+1, c+1);
        }
    }
    _S->FillBandedBlock(order-1, N, [=,&overlapStore](int row, int col) {
        int i, j;
        i = row % N;  
        j = col % N;  

        return overlapStore(i, j);
    });

    //
    if (start_it > 0 && file_exists(_outputFilename)) {
        std::vector<complex> norm((start_it+1)/_computePeriod);
        std::stringstream line;
        
        // first read all the values back that we want to keep
        _txtFile = io::Factory::OpenASCII(_outputFilename, 'r');
        for (int i = 0; i < norm.size(); i++) {
            double real, imag;
            line.str(_txtFile->ReadLine());
            line >> real;
            line >> imag;
            norm[i] = complex(real, imag);
        }

        // now clear the file and write them all back
        _txtFile = io::Factory::OpenASCII(_outputFilename, 'w');
        for (int i = 0; i < norm.size(); i++) {
            line.str("");
            line << std::real(norm[i]) << "\t" << std::imag(norm[i]) << std::endl;
            _txtFile->Write(line.str().c_str());
        }
    } else {
        _txtFile = io::Factory::OpenASCII(_outputFilename, 'w');
    }
}
void NormObservable::Shutdown() {
    _psi = nullptr;
    _psi_temp = nullptr;
    _S = nullptr;
    _txtFile = nullptr;
}
void NormObservable::Compute(int it) {
    complex norm;
    Mult(_S, _psi, _psi_temp);
    Dot(_psi, _psi_temp, norm);
    if (_txtFile) {
        std::stringstream ss;
        ss << std::setprecision(8) << std::scientific;
        ss << std::real(norm) << "\t" << std::imag(norm) << "\n";
        _txtFile->Write(ss.str().c_str());
    } else {
        std::stringstream ss;
        ss << "Norm: " << norm << "\n";
        LOG_INFO(ss.str());
    }
}

void NormObservable::Flush() {
    _txtFile->Flush();
}

std::string NormObservable::GetName() {
    return "norm";
}
Observable::Ptr_t NormObservable::Create(const nlohmann::json& observable) {
    auto norm_obs = new NormObservable();
    if (observable.contains("compute_period")) 
        norm_obs->_computePeriod = observable["compute_period"].get<int>();
    norm_obs->_outputFilename = observable["filename"].get<std::string>();
    return Observable::Ptr_t(norm_obs);
}
bool NormObservable::Validate(const nlohmann::json& observable) {
    if (observable.contains("compute_period") && !observable["compute_period"].is_number()) {
        LOG_CRITICAL("Optional entry \"compute_period\" must be a number.");
        return false;
    }
    if (!observable.contains("filename")) {
        LOG_CRITICAL("\"norm\" observable must contain string entry: filename");
        return false;
    }
    return true;
}