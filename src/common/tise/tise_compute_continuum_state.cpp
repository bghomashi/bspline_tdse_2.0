#include "common/tise/tise.h"


#include "common/file_io/io_factory.h"
#include "common/system_state/system_state.h"
#include "common/utility/logger.h"
#include "common/objects/potential.h"
#include "common/bspline/bspline.h"
#include "objects/potentials/coulomb/coulomb_pot.h"

#include <string>
#include <sstream>
#include <complex>
#include <iomanip>

using namespace tdse;
using namespace maths;

void tise::TISE::ComputeAndOutputContinuumStates(int nmax, Matrix H0, Matrix S, int l) {
    std::stringstream ss;
    auto outputFilename = SystemState::GetEigenStateFilename();
    EigenProblemType pt = SystemState::GetProblemType();
    auto hdf5 = io::Factory::OpenHDF5(outputFilename, 'a');

    _eigensolver->SetProblemType(pt);
    _eigensolver->Solve(H0, S, nmax, _tol);
    auto& values = _eigensolver->GetEigenValues();
    auto& vectors = _eigensolver->GetEigenVectors();



    // find the last negative eigenvalue
    int firstPositiveEigenvalue = 0;
    for (firstPositiveEigenvalue = 0; firstPositiveEigenvalue < values.size(); firstPositiveEigenvalue++)
        if (std::real(values[firstPositiveEigenvalue]) >= 0)
            break;

    // choose a normalization
    if (_continuum_normalization == ContinuumNormalization::MAX) {
        maths::Vector temp = maths::Factory::CreateVector(vectors[0]->Length());

        for (int j = firstPositiveEigenvalue; j < values.size(); j++) {
            auto& v = vectors[j];
            temp->Copy(v);
            temp->Abs();
            double max = temp->Max();
            v->Scale(1./max);
        }
    } else if (_continuum_normalization == ContinuumNormalization::ASYMPTOTICALLY_ONE) {
        // we choose an Rvalue near the edge of the box
        double r = bspline::Basis::GetXmax()*.95;
        std::vector<maths::complex> fc(vectors[0]->Length());
        double Z = 0;
        for (auto& p : _potentials)
            if (p->Name() == CoulombPotential::GetName())
                Z += std::dynamic_pointer_cast<CoulombPotential>(p)->Z();
                

        for (int j = firstPositiveEigenvalue; j < values.size(); j++) {
            double k = std::sqrt(2.*std::real(values[j]));
            auto& v = vectors[j];
            v->CopyTo(fc);
            
            maths::complex phi = bspline::Basis::FunctionEvaluate(r, fc);
            maths::complex dphi = bspline::Basis::FunctionEvaluate(r, fc, 1);

            double A = std::sqrt(
                std::abs(phi)*std::abs(phi) + 
                std::abs(dphi / (k + Z/(k*r))) *
                std::abs(dphi / (k + Z/(k*r)))
            );

            v->Scale(1./A);
        }
    }

    // ----------- output eigen values
    for (int j = firstPositiveEigenvalue; j < values.size(); j++) {
        ss.str("");  
        ss << "energy " << j+1 << " : (" << std::real(values[j]) << " ," << std::imag(values[j]) << ")";
        LOG_INFO(ss.str());
    }

    // store these vectors
    if (l == -1)
        hdf5->WriteAttribute("symmetry", Symmetry::Axial);
    else
        hdf5->WriteAttribute("symmetry", Symmetry::Central);
    hdf5->PushGroup("continuum");
    for (int index = firstPositiveEigenvalue; index < values.size(); index++) {            // n-quantum number
        ss.str("");                                 // clear string stream
        if (l == -1)
            ss << "(" << (index-firstPositiveEigenvalue)+1 << ")";
        else
            ss << "(" << (index-firstPositiveEigenvalue)+1 << ", " << l << ")";

        hdf5->WriteVector(ss.str().c_str(), vectors[index]);
        hdf5->WriteAttribute(vectors[index], "energy", values[index]);
    }
    hdf5->PopGroup();
}