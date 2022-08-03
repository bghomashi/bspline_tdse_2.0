#include "objects/observables/population_obs.h"
#include "common/utility/logger.h"
#include "common/utility/index_manip.h"
#include "common/utility/banded_matrix.h"
#include "common/tdse/simulation.h"
#include "common/file_io/io_factory.h"
#include "common/maths/math_factory.h"
#include "common/maths/math_algebra.h"
#include "common/bspline/bspline.h"
#include "common/system_state/system_state.h"
#include <iostream>
#include <sstream>
#include <iomanip>
#include <string>

using namespace tdse;
using Math = maths::Factory;
using namespace maths;

PopulationObservable::PopulationObservable() {}
void PopulationObservable::Startup(int start_it) {}
void PopulationObservable::Flush() {}
void PopulationObservable::Compute(int it) {/*compute during propagation?*/}



void PopulationObservable::Shutdown() {
    int potSymmetry = Simulation::GetPotentialSymmetry();
    if (potSymmetry == Symmetry::Central)
        PopulationsCentral();
    else if (potSymmetry == Symmetry::Axial)
        PopulationsAxial();
}




std::string PopulationObservable::GetName() {
    return "populations";
}
Observable::Ptr_t PopulationObservable::Create(const nlohmann::json& observable) {
    auto pop_obs = new PopulationObservable();
    // if (observable.contains("compute_period")) 
    //     norm_obs->_computePeriod = observable["compute_period"].get<int>();
    pop_obs->_outputFilename = observable["filename"].get<std::string>();
    return Observable::Ptr_t(pop_obs);
}
bool PopulationObservable::Validate(const nlohmann::json& observable) {
    // compute during propagation?
    // if (observable.contains("compute_period") && !observable["compute_period"].is_number()) {
    //     LOG_CRITICAL("Optional entry \"compute_period\" must be a number.");
    //     return false;
    // }
    if (!observable.contains("filename")) {
        LOG_CRITICAL("\"populations\" observable must contain string entry: filename");
        return false;
    }
    return true;
}


void PopulationObservable::PopulationsCentral() {
    complex pop;
    std::stringstream ss;
    io::ASCII txtFile;
    // ---------- grab a bunch of state variables ------
    std::string eigenStateFilename = SystemState::GetEigenStateFilename();
    int eigenStateLmax = SystemState::GetEigenStateBoundLmax();
    int eigenStateBoundNmax = SystemState::GetEigenStateBoundNmax();
    int lmax = SystemState::GetBasisLmax();
    int order = bspline::Basis::GetOrder();
    int N = bspline::Basis::GetNumBSplines();
    auto& Ms = Simulation::GetMs();
    auto& mRows = Simulation::GetMRows();
    auto& psi = Simulation::GetPsi();
    auto hdf5 = io::Factory::OpenHDF5(eigenStateFilename, 'r');
    Vector psiTemp = Math::CreateVector(N);
    Vector eigenState = maths::Factory::CreateVector(N);
    Matrix S = Math::CreateMatrix(N, N, 2*order-1);
    ss << std::setprecision(8) << std::scientific;

    // --------- build overlap -------------
    BuildOverlap(S);

    // ----------- open text file
    if ((txtFile = io::Factory::OpenASCII(_outputFilename, 'w')) == nullptr) {
        LOG_INFO("Failed to open file.");
        return;
    }
    ss << "n\tl\tm";
    ss << std::setw(20) << std::right << "pop";
    ss << std::setw(20) << std::right << "real(Amp)";
    ss << std::setw(20) << std::right << "imag(Amp)";
    ss << "\n";
    txtFile->Write(ss.str());
    hdf5->PushGroup("bound");

    for (auto m : Ms) {                             // for each M state in Psi
        if (std::abs(m) > eigenStateLmax)           // we do not have this state to project on to
            continue;                               // so skip

        for (int l = std::abs(m); l <= std::min(lmax, eigenStateLmax); l++) { // for each L in this M
            // get a subvector
            int start = RowFrom(m, Ms, mRows) + (l-std::abs(m))*N;
            Vector lm_block = psi->GetSubVector(start, start + N);
            
            for (int n = l+1; n <= eigenStateBoundNmax; n++) {           // for each N in this LM state
                // load vector from state file
                ss.str("");                                         // clear string stream
                ss << "(" << n << ", " << l << ")";                 // name of state
                hdf5->ReadVector(ss.str().c_str(), eigenState);     

                // ------------ project ----------------
                Mult(S, eigenState, psiTemp);
                Dot(lm_block, psiTemp, pop);
                
                // ----------- output ----------------
                ss.str("");
                ss << n << "\t" << l << "\t" << m << "\t";
                ss << std::abs(pop)*std::abs(pop) << "\t" 
                << std::real(pop) << "\t" 
                << std::imag(pop) << "\n";

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
void PopulationObservable::PopulationsAxial() {
    complex pop;
    std::stringstream ss;
    io::ASCII txtFile;
    std::string eigenStateFilename = SystemState::GetEigenStateFilename();
    int lmax = SystemState::GetBasisLmax();
    int eigenStateBoundNmax = SystemState::GetEigenStateBoundNmax();
    int order = bspline::Basis::GetOrder();
    int N = bspline::Basis::GetNumBSplines();
    auto& Ms = Simulation::GetMs();
    auto& mRows = Simulation::GetMRows();
    auto& psi = Simulation::GetPsi();
    auto hdf5 = io::Factory::OpenHDF5(eigenStateFilename, 'r');
    Vector psiTemp = Math::CreateVector(N*(lmax+1));
    Vector eigenState = maths::Factory::CreateVector(N*(lmax+1));
    Matrix S = Math::CreateMatrix(N*(lmax+1), N*(lmax+1), 2*order-1);

    ss << std::setprecision(8) << std::scientific;

    // --------- build overlap -------------
    BuildOverlap(S);

    // ----------- open text file
    if ((txtFile = io::Factory::OpenASCII(_outputFilename, 'w')) == nullptr) {
        LOG_INFO("Failed to open file.");
        return;
    }

    ss << "n\tm";
    ss << std::setw(20) << std::right << "pop";
    ss << std::setw(20) << std::right << "real(Amp)";
    ss << std::setw(20) << std::right << "imag(Amp)";
    ss << "\n";
    txtFile->Write(ss.str());
    hdf5->PushGroup("bound");

    for (auto m : Ms) {                             // for each M state in Psi
        // get a subvector
        int start = RowFrom(m, Ms, mRows);
        Vector mBlock = psi->GetSubVector(start, start + N*(lmax+1));
        
        for (int n = 1; n <= eigenStateBoundNmax; n++) {           // for each N in this LM state
            // load vector from state file
            ss.str("");                                         // clear string stream
            ss << "(" << n << ")";                              // name of state
            hdf5->ReadVector(ss.str().c_str(), eigenState);     

            // ------------ project ----------------
            Mult(S, eigenState, psiTemp);
            Dot(mBlock, psiTemp, pop);
            
            // ----------- output ----------------
            ss.str("");
            ss << n << "\t" << m << "\t";
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

void PopulationObservable::BuildOverlap(Matrix S) {
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