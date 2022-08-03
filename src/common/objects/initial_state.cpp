#include "common/objects/initial_state.h"
#include "common/tdse/simulation.h"
#include "common/maths/math_factory.h"
#include "common/maths/vector.h"
#include "common/maths/math_algebra.h"
#include "common/file_io/io_factory.h"
#include "common/utility/file_exists.h"
#include "common/utility/logger.h"
#include "common/utility/index_manip.h"
#include "common/bspline/bspline.h"
#include "common/system_state/system_state.h"



#include "petsc/maths/petsc_common.h"

using namespace std::complex_literals;
using namespace tdse;


void InitialState::SetFilename(const std::string& filename) {
    _eigenStateFilename = filename;
}
void InitialState::BuildInitialState(maths::Vector psi) {
    if (!file_exists(_eigenStateFilename)) {
        LOG_CRITICAL("eigenstate file does not exists: " + _eigenStateFilename);
        exit(-1);
    }

    int lmax = SystemState::GetBasisLmax();
    int dof = Simulation::GetDOF();                       // degrees of freedom
    auto Ms = Simulation::GetMs();
    auto MRows = Simulation::GetMRows();
    int N = bspline::Basis::GetNumBSplines();               // number of bsplines
    int potSymmetry;
    double norm = 0;

    // compute norm from amplitudes
    for (auto& state : _initialState)
        norm += state.amplitude*state.amplitude;

    // ------------- open file ----------------
    auto hdf5 = io::Factory::OpenHDF5(_eigenStateFilename, 'r');
    hdf5->ReadAttribute("symmetry", &potSymmetry);
    hdf5->PushGroup("bound");

    if (potSymmetry == Symmetry::Central) {
        maths::Vector temp = maths::Factory::CreateVector(N);   // vector from hdf5 file
        std::vector<maths::Vector> vecs(dof/N);                 // one vector for each lm - these will all be concatenated
        for (auto& v : vecs) {
            v = maths::Factory::CreateVector(N); 
            v->Zero();
        }


        // merge the initial eigenstates into 'vecs'
        std::stringstream name_ss;
        for (auto& state : _initialState) {
            name_ss.str("");                                            // clear string stream
            name_ss << "(" << state.n << ", " << state.l << ")";        // name of state

            int mBlock = RowFrom(state.m, Ms, MRows)/N;
            int lmBlock = mBlock + (state.l-std::abs(state.m));
            hdf5->ReadVector(name_ss.str().c_str(), temp);              // read in the state
            temp->Scale(state.amplitude*std::exp(1.i*state.phase));     // scale by amplitude and phase
            maths::AXPY(vecs[lmBlock], 1., temp);                        // sum with other similar lm's
        }
        psi->Concatenate(vecs);                                         // append all the initial vectors together

    } else if (potSymmetry == Symmetry::Axial) {
        maths::Vector temp = maths::Factory::CreateVector(N*(lmax+1));           // vector from hdf5 file
        std::vector<maths::Vector> vecs(dof/(N*(lmax+1)));                 // one vector for each m
        for (auto& v : vecs) {
            v = maths::Factory::CreateVector(N*(lmax+1)); 
            v->Zero();
        }

        std::stringstream name_ss;
        for (auto& state : _initialState) {
            name_ss.str("");                                            // clear string stream
            name_ss << "(" << state.n << ")";                           // name of state

            int mBlock = RowFrom(state.m, Ms, MRows)/(N*(lmax+1));          
            hdf5->ReadVector(name_ss.str().c_str(), temp);              // read in the state
            temp->Scale(state.amplitude*std::exp(1.i*state.phase));     // scale by amplitude and phase
            maths::AXPY(vecs[mBlock], 1., temp);                        // sum with other similar lm's
        }
        psi->Concatenate(vecs);                                         // append all the initial vectors together
    }
    hdf5->PopGroup();

    // concatenate into "psi"
    psi->Scale(1./std::sqrt(norm));                     // and normalize
}
void InitialState::AddInitialState(int n, int l, int m, double phase,  double amplitude) {
    _initialState.emplace_back(StateDescriptor{n, l, m, phase, amplitude});
}


bool InitialState::Validate(const nlohmann::json& input) {
    for (auto& state : input) {
        if (!(state.contains("n") && state["n"].is_number())) {
            LOG_CRITICAL("initial state must contain number entry: n");
            return false;
        }
        if (!(state.contains("l") && state["l"].is_number())) {
            LOG_CRITICAL("initial state must contain number entry: l");
            return false;
        }
        if (!(state.contains("m") && state["m"].is_number())) {
            LOG_CRITICAL("initial state must contain number entry: m");
            return false;
        }
        if (state.contains("phase") && !state["phase"].is_number()) {
            LOG_CRITICAL("Optional entry \"phase\" must be a number.");
            return false;
        }
        if (state.contains("amplitude") && !state["amplitude"].is_number()) {
            LOG_CRITICAL("Optional entry \"amplitude\" must be a number.");
            return false;
        }
    }
    return true;
}
void InitialState::Load(const nlohmann::json& input) {
    for (auto& state : input) {
        int n = state["n"];
        int l = state["l"];
        int m = state["m"];
        double phase = 0.;
        double amplitude = 1.;

        if (state.contains("phase")) phase = state["phase"];
        if (state.contains("amplitude")) amplitude = state["amplitude"];
        
        AddInitialState(n, l, m, phase, amplitude);
    }
}