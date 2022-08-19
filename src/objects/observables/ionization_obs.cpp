#include "objects/observables/ionization_obs.h"
#include "common/utility/logger.h"
#include "common/utility/index_manip.h"
#include "common/utility/banded_matrix.h"
#include "common/utility/gsl_fit_sin.h"
#include "common/tdse/simulation.h"
#include "common/file_io/io_factory.h"
#include "common/maths/math_factory.h"
#include "common/maths/math_algebra.h"
#include "common/bspline/bspline.h"
#include "common/system_state/system_state.h"
#include "objects/potentials/coulomb/coulomb_pot.h"
#include <iostream>
#include <sstream>
#include <iomanip>
#include <string>

using namespace tdse;
using Math = maths::Factory;
using namespace maths;
using namespace std::complex_literals;

IonizationObservable::IonizationObservable() {}
void IonizationObservable::Startup(int start_it) {}
void IonizationObservable::Flush() {}
void IonizationObservable::Compute(int it) {}



void IonizationObservable::Shutdown() {
    int potSymmetry = Simulation::GetPotentialSymmetry();
    if (potSymmetry == Symmetry::Central)
        IonizationCentral();
    else if (potSymmetry == Symmetry::Axial)
        IonizationAxial();
}




std::string IonizationObservable::GetName() {
    return "ionization";
}
Observable::Ptr_t IonizationObservable::Create(const nlohmann::json& observable) {
    auto ion_obs = new IonizationObservable();
    ion_obs->_outputFilename = observable["filename"].get<std::string>();
    return Observable::Ptr_t(ion_obs);
}
bool IonizationObservable::Validate(const nlohmann::json& observable) {
    if (!observable.contains("filename")) {
        LOG_CRITICAL("\"ionization\" observable must contain string entry: filename");
        return false;
    }
    return true;
}


void IonizationObservable::IonizationCentral() {
    complex pop;
    std::stringstream ss;
    io::ASCII txtFile;
    // ---------- grab a bunch of state variables ------
    std::string eigenStateFilename = SystemState::GetEigenStateFilename();
    int eigenStateContinuumLmax = SystemState::GetEigenStateContinuumLmax();
    int eigenStateContinuumNmax = SystemState::GetEigenStateContinuumNmax();
    int lmax = SystemState::GetBasisLmax();
    int order = bspline::Basis::GetOrder();
    int N = bspline::Basis::GetNumBSplines();
    auto& potentials = Simulation::GetPotentials();
    auto& Ms = Simulation::GetMs();
    auto& mRows = Simulation::GetMRows();
    auto& psi = Simulation::GetPsi();
    auto hdf5 = io::Factory::OpenHDF5(eigenStateFilename, 'r');
    Vector psiTemp = Math::CreateVector(N);
    Vector eigenState = maths::Factory::CreateVector(N);
    Matrix S = Math::CreateMatrix(N, N, 2*order-1);
    double energy;

    // for fitting phaseshift
    double Z = 0;
    for (auto& p : potentials)
        if (p->Name() == CoulombPotential::GetName())
            Z += std::dynamic_pointer_cast<CoulombPotential>(p)->Z();

    double rmin = bspline::Basis::GetXmax()*.9;
    double rmax = bspline::Basis::GetXmax()*.95;
    int nr = 100;
    double dr = (rmax - rmin) / (nr - 1.);
    std::vector<double> out(2);
    std::vector<double> rs(nr), ys(nr);
    std::vector<maths::complex> state(N), ys_complex(nr);
    for (int i = 0; i < nr; i++)
        rs[i] = rmin + i*dr;
    //

    ss << std::setprecision(8) << std::scientific;

    // --------- build overlap -------------
    BuildOverlap(S);

    // ----------- open text file
    if ((txtFile = io::Factory::OpenASCII(_outputFilename, 'w')) == nullptr) {
        LOG_CRITICAL("Failed to open file.");
        return;
    }
    ss << std::setw(14) << std::right << "E\t";
    ss << "l\tm";
    ss << std::setw(15) << std::right << "pop";
    ss << std::setw(15) << std::right << "real(Amp)";
    ss << std::setw(15) << std::right << "imag(Amp)";
    ss << std::setw(15) << std::right << "phaseshift";
    ss << "\n";
    txtFile->Write(ss.str());
    hdf5->PushGroup("continuum");

    for (auto m : Ms) {                             // for each M state in Psi
        if (std::abs(m) > eigenStateContinuumLmax)           // we do not have this state to project on to
            continue;                               // so skip

        for (int l = std::abs(m); l <= std::min(lmax, eigenStateContinuumLmax); l++) { // for each L in this M
            // get a subvector
            int start = RowFrom(m, Ms, mRows) + (l-std::abs(m))*N;
            Vector lm_block = psi->GetSubVector(start, start + N);
            
            for (int n = 1; n <= eigenStateContinuumNmax; n++) {           // for each N in this LM state
                // load vector from state file
                ss.str("");                                         // clear string stream
                ss << "(" << n << ", " << l << ")";                 // name of state
                if (!hdf5->HasVector(ss.str().c_str()))
                    continue;
                hdf5->ReadVector(ss.str().c_str(), eigenState);     
                hdf5->ReadAttribute(eigenState, "energy-real", &energy);

                // ---------- get phase --------------
                {
                    
                    eigenState->CopyTo(state);

                    // here another option is to abs-square state
                    ys_complex = bspline::Basis::FunctionEvaluate(rs, state);
                    for (int i = 0; i < ys.size(); i++)
                        ys[i] = std::real(ys_complex[i]);
                    
                    if (!fit_sin({1.0, 0.0}, rs, ys, {sqrt(2.*energy), (double)l, Z}, out))
                        LOG_CRITICAL("failed to get eigenstate phaseshift");

                    // if (out[0] < 0) {
                    //     out[1] += maths::Pi;
                    //     out[0] *= -1;
                    // }
                    // out[0] is the amplitude of the asymptotic eigenstate ;
                    // out[1] is the phaseshift

                }

                // ------------ project ----------------
                Mult(S, eigenState, psiTemp);
                Dot(lm_block, psiTemp, pop);

                pop *= std::exp(1.i*out[1]) / out[0];

                // ----------- output ----------------
                ss.str("");
                ss << energy << "\t" << l << "\t" << m << "\t";
                ss << std::abs(pop)*std::abs(pop) << "\t" 
                << std::real(pop) << "\t"
                << std::imag(pop) << "\t"
                << out[1] << "\n";

                txtFile->Write(ss.str().c_str());
            }
            psi->RestoreSubVector(lm_block);

            txtFile->Write("\n");
        }
    }

    hdf5->PopGroup();

    hdf5 = nullptr;
    psiTemp = nullptr;
    eigenState = nullptr;
    S = nullptr;
}
void IonizationObservable::IonizationAxial() {
    complex pop;
    std::stringstream ss;
    io::ASCII txtFile;
    std::string eigenStateFilename = SystemState::GetEigenStateFilename();
    int lmax = SystemState::GetBasisLmax();
    int eigenStateContinuumNmax = SystemState::GetEigenStateContinuumNmax();
    int order = bspline::Basis::GetOrder();
    int N = bspline::Basis::GetNumBSplines();
    auto& Ms = Simulation::GetMs();
    auto& mRows = Simulation::GetMRows();
    auto& psi = Simulation::GetPsi();
    auto hdf5 = io::Factory::OpenHDF5(eigenStateFilename, 'r');
    Vector psiTemp = Math::CreateVector(N*(lmax+1));
    Vector eigenState = maths::Factory::CreateVector(N*(lmax+1));
    double energy;
    Matrix S = Math::CreateMatrix(N*(lmax+1), N*(lmax+1), 2*order-1);

    ss << std::setprecision(8) << std::scientific;

    // --------- build overlap -------------
    BuildOverlap(S);

    // ----------- open text file
    if ((txtFile = io::Factory::OpenASCII(_outputFilename, 'w')) == nullptr) {
        LOG_INFO("Failed to open file.");
        return;
    }

    ss << std::setw(14) << std::right << "E\t";
    ss << "m";
    ss << std::setw(20) << std::right << "pop";
    ss << std::setw(20) << std::right << "real(Amp)";
    ss << std::setw(20) << std::right << "imag(Amp)";
    ss << "\n";
    txtFile->Write(ss.str());
    hdf5->PushGroup("continuum");

    for (auto m : Ms) {                             // for each M state in Psi
        // get a subvector
        int start = RowFrom(m, Ms, mRows);
        Vector mBlock = psi->GetSubVector(start, start + N*(lmax+1));
        
        for (int n = 1; n <= eigenStateContinuumNmax; n++) {           // for each N in this LM state
            // load vector from state file
            ss.str("");                                         // clear string stream
            ss << "(" << n << ")";                              // name of state
            if (!hdf5->HasVector(ss.str().c_str()))
                break;
            hdf5->ReadVector(ss.str().c_str(), eigenState); 
            hdf5->ReadAttribute(eigenState, "energy-real", &energy);

            // ------------ project ----------------
            Mult(S, eigenState, psiTemp);
            Dot(mBlock, psiTemp, pop);
            
            // ----------- output ----------------
            ss.str("");
            ss << "(" << energy << "\t" << m << "): ";
            ss << std::abs(pop)*std::abs(pop) << "\t" 
            << std::real(pop) << "\t" 
            << std::imag(pop) << "\n";

            txtFile->Write(ss.str().c_str());
        }

        psi->RestoreSubVector(mBlock);
        txtFile->Write("\n");
    }

    hdf5->PopGroup();

    hdf5 = nullptr;
    psiTemp = nullptr;
    eigenState = nullptr;
    S = nullptr;
}

void IonizationObservable::BuildOverlap(Matrix S) {
    int N = bspline::Basis::GetNumBSplines();
    int order = bspline::Basis::GetNumBSplines();
    BandedMatrix overlapStore(N, 2*order-1);

    for (int i = 0; i < N; i++) {
        for (int j = i; j < std::min(N, i+order); j++) {
            overlapStore(j,i) = overlapStore(i,j) = bspline::Basis::Integrate(i+1, j+1);
        }
    }
    S->FillBandedBlock(order-1, N, [=,&overlapStore](int row, int col) {
        int i, j;
        i = row % N;  
        j = col % N;  

        return overlapStore(i, j);
    });
}