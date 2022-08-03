#include "objects/observables/dipole_obs.h"
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include "common/tdse/simulation.h"
#include "common/utility/banded_matrix.h"
#include "common/utility/logger.h"
#include "common/utility/file_exists.h"
#include "common/maths/math_factory.h"
#include "common/maths/math_algebra.h"
#include "common/file_io/io_factory.h"
#include "common/utility/index_manip.h"
#include "common/utility/spherical_harmonics.h"
#include "common/bspline/bspline.h"
#include "common/system_state/system_state.h"

using namespace tdse;
using Math = maths::Factory;
using namespace maths;

DipoleObservable::DipoleObservable() {
}
void DipoleObservable::Startup(int start_it) {
    // Build GradPotential Matrix
    int Nmax = bspline::Basis::GetNumBSplines();
    int order = bspline::Basis::GetOrder();
    int Lmax = SystemState::GetBasisLmax();
    int dof = Simulation::GetDOF();
    auto& Ms = Simulation::GetMs();
    auto& mRows = Simulation::GetMRows();
    auto polarization = Simulation::GetPolarization();
    
    _psi = Simulation::GetPsi();
    _psiTemp = Math::CreateVector(_psi->Length());

    if (polarization[X]) {
        Matrix temp = Math::CreateMatrix(dof, dof, 8*order-4);
        _x[X] = Math::CreateMatrix(dof, dof, 8*order-4);

        // Fill
        FillMatrixX(_x[X], Nmax, Lmax, Ms, mRows);
    }
    if (polarization[Y]) {
        Matrix temp = Math::CreateMatrix(dof, dof, 8*order-4);
        _x[Y] = Math::CreateMatrix(dof, dof, 8*order-4);

        // Fill
        FillMatrixY(_x[Y], Nmax, Lmax, Ms, mRows);
    }
    if (polarization[Z]) {
        Matrix temp = Math::CreateMatrix(dof, dof, 4*order-2);
        _x[Z] = Math::CreateMatrix(dof, dof, 4*order-2);

        // Fill
        FillMatrixZ(_x[Z], Nmax, Lmax, Ms, mRows);
    }
    
    // here we want to clear any extra entrees
    if (start_it > 0 && file_exists(_outputFilename)) {
        std::vector<double> t((start_it+1)/_computePeriod), x((start_it+1)/_computePeriod), y((start_it+1)/_computePeriod), z((start_it+1)/_computePeriod);
        std::stringstream line;
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

void DipoleObservable::Shutdown() {
    _psi = nullptr;
    _psiTemp = nullptr;
    _x[X] = nullptr;
    _x[Y] = nullptr;
    _x[Z] = nullptr;
    _txtFile = nullptr;
}


void DipoleObservable::Compute(int it) {
    double t = Simulation::GetTime();
    complex dipole[DimIndex::NUM] = {0., 0., 0.};
    std::stringstream ss;

    if (_x[X]) {
        Mult(_x[X], _psi, _psiTemp);
        Dot(_psi, _psiTemp, dipole[X]);
    }
    if (_x[Y]) {
        Mult(_x[Y], _psi, _psiTemp);
        Dot(_psi, _psiTemp, dipole[Y]);
    }
    if (_x[Z]) {
        Mult(_x[Z], _psi, _psiTemp);
        Dot(_psi, _psiTemp, dipole[Z]);
    }
    
    ss << std::setprecision(8) << std::scientific;
    ss  << t << "\t" 
        << std::real(dipole[X]) << "\t"
        << std::real(dipole[Y]) << "\t"
        << std::real(dipole[Z]) << std::endl;
    _txtFile->Write(ss.str().c_str());
}

void DipoleObservable::Flush() {
    _txtFile->Flush();
}

int DipoleObservable::MemoryAlloced() const {
    int order = bspline::Basis::GetOrder();
    int dof = Simulation::GetDOF();

    int mem = dof;
    if (_x[X])
        mem += dof*(8*order-4);
    if (_x[Y])
        mem += dof*(8*order-4);
    if (_x[Z])
        mem += dof*(4*order-2);
    
    return mem;
}

std::string DipoleObservable::GetName() {
    return "dipole";
}
Observable::Ptr_t DipoleObservable::Create(const nlohmann::json& observable) {
    auto dip_acc_obs = new DipoleObservable();
    if (observable.contains("compute_period")) 
        dip_acc_obs->_computePeriod = observable["compute_period"].get<int>();
    dip_acc_obs->_outputFilename = observable["filename"].get<std::string>();
    return Observable::Ptr_t(dip_acc_obs);
}
bool DipoleObservable::Validate(const nlohmann::json& observable) {
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
void DipoleObservable::BuildR(int N, int order, BandedMatrix& invR) {
    for (int i = 0; i < N; i++) {
        for (int j = i; j < std::min(N,i+order); j++) {
            invR(i, j) = invR(j, i) = bspline::Basis::Integrate(i+1, j+1, [] (complex r) -> complex {
                return r;
            });
        }
    }
}
void DipoleObservable::FillBlock( int l1, int m1, int l2, int m2,
                int lmax, int mmax, int N, int order,
                const std::vector<int>& Ms,
                const std::vector<int>& mRows,
                const BandedMatrix& r,
                std::function<complex(int, int, int, int)> YlmYlm,
                Matrix m) {
    if (l2 <= lmax && l2 >= std::abs(m2) &&
        l2 >= 0 && l2 >= std::abs(m2) &&
        m2 <= mmax && m2 >= -mmax) {
        int m1Block = RowFrom(m1, Ms, mRows)/N;
        int m2Block = RowFrom(m2, Ms, mRows)/N;
        int blockRow = m1Block + (l1-std::abs(m1));
        int blockCol = m2Block + (l2-std::abs(m2));

        m->FillBandedBlock(order-1, N, blockRow, blockCol, 
        [=,&r](int row, int col) {
            int i = row % N, j = col % N;
            return r(i, j)*YlmYlm(l1,m1,l2,m2);
        });
    }
}
void DipoleObservable::FillMatrixX(Matrix m, int N, int lmax, const std::vector<int>& Ms, const std::vector<int>& mRows) {
    int order = bspline::Basis::GetOrder();
    // if we get here this is a full 3D calculations
    // so Ms = [-mmax,mmax]
    int mmax = Ms.back();

    // the same for each L-M block so cache it
    BandedMatrix R(N,2*order-1);
    BuildR(N, order, R);

    for (int m1 : Ms) {
        for (int l1 = std::abs(m1); l1 <= lmax; l1++) {
            FillBlock(l1, m1, l1+1, m1+1,
                        lmax, mmax, N, order, 
                        Ms, mRows, R, YlmXYlm, m);
            FillBlock(l1, m1, l1-1, m1+1,
                        lmax, mmax, N, order, 
                        Ms, mRows, R, YlmXYlm, m);
            FillBlock(l1, m1, l1+1, m1-1,
                        lmax, mmax, N, order, 
                        Ms, mRows, R, YlmXYlm, m);
            FillBlock(l1, m1, l1-1, m1-1,
                        lmax, mmax, N, order, 
                        Ms, mRows, R, YlmXYlm, m);
        }
    }

    m->AssembleBegin();
    m->AssembleEnd();
}

void DipoleObservable::FillMatrixY(Matrix m, int N, int lmax, const std::vector<int>& Ms, const std::vector<int>& mRows) {
    int order = bspline::Basis::GetOrder();
    int mmax = Ms.back();

    // the same for each L-M block so cache it
    BandedMatrix R(N,2*order-1);
    BuildR(N, order, R);

    // for each m-block (block rows)
    for (int m1 : Ms) {
        for (int l1 = std::abs(m1); l1 <= lmax; l1++) {
            FillBlock(l1, m1, l1+1, m1+1,
                        lmax, mmax, N, order, 
                        Ms, mRows, R, YlmYYlm, m);
            FillBlock(l1, m1, l1-1, m1+1,
                        lmax, mmax, N, order, 
                        Ms, mRows, R, YlmYYlm, m);
            FillBlock(l1, m1, l1+1, m1-1,
                        lmax, mmax, N, order, 
                        Ms, mRows, R, YlmYYlm, m);
            FillBlock(l1, m1, l1-1, m1-1,
                        lmax, mmax, N, order, 
                        Ms, mRows, R, YlmYYlm, m);
        }
    }

    m->AssembleBegin();
    m->AssembleEnd();
}

void DipoleObservable::FillMatrixZ(Matrix m, int N, int lmax, const std::vector<int>& Ms, const std::vector<int>& mRows) {
    int order = bspline::Basis::GetOrder();
    int mmax = Ms.back();

    // the same for each L-M block so cache it
    BandedMatrix R(N, 2*order-1);
    BuildR(N, order, R);
    for (int m1 : Ms) {
        for (int l1 = std::abs(m1); l1 <= lmax; l1++) { 
            FillBlock(l1, m1, l1+1, m1,
                        lmax, mmax, N, order, 
                        Ms, mRows, R, YlmZYlm, m);
            FillBlock(l1, m1, l1-1, m1,
                        lmax, mmax, N, order, 
                        Ms, mRows, R, YlmZYlm, m);
        }
    }
    m->AssembleBegin();
    m->AssembleEnd();
}