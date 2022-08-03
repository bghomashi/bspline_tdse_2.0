#include "objects/observables/dipole_acc_obs.h"
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include "common/tdse/simulation.h"
#include "common/utility/logger.h"
#include "common/utility/file_exists.h"
#include "common/maths/math_factory.h"
#include "common/maths/math_algebra.h"
#include "common/file_io/io_factory.h"
#include "common/bspline/bspline.h"
#include "common/system_state/system_state.h"

using namespace tdse;
using Math = maths::Factory;
using namespace maths;

DipoleAccObservable::DipoleAccObservable() {
}
void DipoleAccObservable::Startup(int start_it) {
    // Build GradPotential Matrix
    int Nmax = bspline::Basis::GetNumBSplines();
    int order = bspline::Basis::GetOrder();
    int Lmax = SystemState::GetBasisLmax();
    int dof = Simulation::GetDOF();
    auto& Ms = Simulation::GetMs();
    auto& mRows = Simulation::GetMRows();
    auto& potentials = Simulation::GetPotentials();
    auto polarization = Simulation::GetPolarization();
    
    _psi = Simulation::GetPsi();
    _psiTemp = Math::CreateVector(_psi->Length());

    if (polarization[X]) {
        Matrix temp = Math::CreateMatrix(dof, dof, 8*order-4);
        _gradPot[X] = Math::CreateMatrix(dof, dof, 8*order-4);

        // Fill
        potentials[0]->FillMatrixGradX(_gradPot[X], Nmax, Lmax, Ms, mRows);
        _gradPot[X]->Scale(-1.);
        for (int i = 1; i < potentials.size(); i++) {
            const auto& p = potentials[i]; 
            p->FillMatrixGradX(temp, Nmax, Lmax, Ms, mRows);    // get potential matrix
            AXPY(_gradPot[X], -1., temp);                   // add potential to H0
        }
    }
    if (polarization[Y]) {
        Matrix temp = Math::CreateMatrix(dof, dof, 8*order-4);
        _gradPot[Y] = Math::CreateMatrix(dof, dof, 8*order-4);

        // Fill
        potentials[0]->FillMatrixGradY(_gradPot[Y], Nmax, Lmax, Ms, mRows);
        _gradPot[Y]->Scale(-1.);
        for (int i = 1; i < potentials.size(); i++) {
            const auto& p = potentials[i]; 
            p->FillMatrixGradY(temp, Nmax, Lmax, Ms, mRows);    // get potential matrix
            AXPY(_gradPot[Y], -1., temp);                   // add potential to H0
        }
    }
    if (polarization[Z]) {
        Matrix temp = Math::CreateMatrix(dof, dof, 4*order-2);
        _gradPot[Z] = Math::CreateMatrix(dof, dof, 4*order-2);

        // Fill
        potentials[0]->FillMatrixGradZ(_gradPot[Z], Nmax, Lmax, Ms, mRows);
        _gradPot[Z]->Scale(-1.);
        for (int i = 1; i < potentials.size(); i++) {
            const auto& p = potentials[i]; 
            p->FillMatrixGradZ(temp, Nmax, Lmax, Ms, mRows);    // get potential matrix
            AXPY(_gradPot[Z], -1., temp);                   // add potential to H0
        }
    }
    
    // here we want to clear any extra entrees
    if (start_it > 0 && file_exists(_outputFilename)) {
        int numKeeping = start_it/_computePeriod + 1;
        std::vector<double> t(numKeeping), x(numKeeping), y(numKeeping), z(numKeeping);
        std::stringstream line;
        line << std::setprecision(8) << std::scientific;
        // first read all the values back that we want to keep
        _txtFile = io::Factory::OpenASCII(_outputFilename, 'r');
        for (int i = 0; i < t.size(); i++) {
            line.str(_txtFile->ReadLine());
            line >> t[i] >> x[i] >> y[i] >> z[i];
        }

        // now clear the file and write them all back
        _txtFile = io::Factory::OpenASCII(_outputFilename, 'w');
        for (int i = 0; i < t.size(); i++) {
            line.str("");
            line << t[i] << "\t" 
                << x[i] << "\t"
                << y[i] << "\t"
                << z[i] << std::endl;
            _txtFile->Write(line.str().c_str());
        }
    } else {
        _txtFile = io::Factory::OpenASCII(_outputFilename, 'w');
    }
}

void DipoleAccObservable::Shutdown() {
    _psi = nullptr;
    _psiTemp = nullptr;
    _gradPot[X] = nullptr;
    _gradPot[Y] = nullptr;
    _gradPot[Z] = nullptr;
    _txtFile = nullptr;
}


void DipoleAccObservable::Compute(int it) {
    double t = Simulation::GetTime();
    complex dipole[DimIndex::NUM] = {0., 0., 0.};
    std::stringstream ss;

    if (_gradPot[X]) {
        Mult(_gradPot[X], _psi, _psiTemp);
        Dot(_psi, _psiTemp, dipole[X]);
    }
    if (_gradPot[Y]) {
        Mult(_gradPot[Y], _psi, _psiTemp);
        Dot(_psi, _psiTemp, dipole[Y]);
    }
    if (_gradPot[Z]) {
        Mult(_gradPot[Z], _psi, _psiTemp);
        Dot(_psi, _psiTemp, dipole[Z]);
    }
    
    ss << std::setprecision(8) << std::scientific;
    ss  << t << "\t" 
        << std::real(dipole[X]) << "\t"
        << std::real(dipole[Y]) << "\t"
        << std::real(dipole[Z]) << std::endl;
    _txtFile->Write(ss.str().c_str());
}

void DipoleAccObservable::Flush() {
    _txtFile->Flush();
}

int DipoleAccObservable::MemoryAlloced() const {
    int order = bspline::Basis::GetOrder();
    int dof = Simulation::GetDOF();

    int mem = dof;
    if (_gradPot[X])
        mem += dof*(8*order-4);
    if (_gradPot[Y])
        mem += dof*(8*order-4);
    if (_gradPot[Z])
        mem += dof*(4*order-2);
    
    return mem;
}

std::string DipoleAccObservable::GetName() {
    return "dipole_acc";
}
Observable::Ptr_t DipoleAccObservable::Create(const nlohmann::json& observable) {
    auto dip_acc_obs = new DipoleAccObservable();
    if (observable.contains("compute_period")) 
        dip_acc_obs->_computePeriod = observable["compute_period"].get<int>();
    dip_acc_obs->_outputFilename = observable["filename"].get<std::string>();
    return Observable::Ptr_t(dip_acc_obs);
}
bool DipoleAccObservable::Validate(const nlohmann::json& observable) {
    // an output file is required
    if (observable.contains("compute_period") && !observable["compute_period"].is_number()) {
        LOG_CRITICAL("Optional entry \"compute_period\" must be a number.");
        return false;
    }
    // an output file is required
    if (!observable.contains("filename")) {
        LOG_CRITICAL("\"dipole_acc\" observable must contain string entry: filename");
        return false;
    }
    return true;
}
