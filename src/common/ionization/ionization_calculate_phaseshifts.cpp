#include "common/ionization/ionization.h"
#include "common/maths/math_common.h"
#include "common/maths/math_factory.h"
#include "common/maths/vector.h"
#include "common/file_io/io_factory.h"
#include "common/bspline/bspline.h"
#include <sstream>
#include <string>

using namespace std::complex_literals;

double phase_shift(double k, int l, double Z, double r, maths::complex phi, maths::complex dphi);

void Ionization::CalculatePhaseshifts() {
    std::stringstream ss;
    double energy, Rmax = bspline::Basis::GetXmax()*.9;
    double Z = 2;
    int N = bspline::Basis::GetNumBSplines();
    auto hdf5 = io::Factory::OpenHDF5(_eigenstate_filename, 'r');
    maths::Vector eigenState = maths::Factory::CreateVector(N);
    std::vector<maths::complex> array(N);
    

    hdf5->PushGroup("continuum");

    // for each continuum eigenstate
    for (int l = 0; l <= lmax; l++) {
        for (int n = 1; n <= continuum_nmax; n++) {           // for each N in this LM state
            ss.str("");                                         // clear string stream
            ss << "(" << n << ", " << l << ")";                 // name of state
            
            // read it into a vector
            if (!hdf5->HasVector(ss.str().c_str()))
                continue;
            hdf5->ReadVector(ss.str().c_str(), eigenState);     
            hdf5->ReadAttribute(eigenState, "energy-real", &energy);
            
            // convert to std::vector
            eigenState->CopyTo(array);
            maths::complex phi, dphi, max;
            // evaluate at Rmax
            phi = bspline::Basis::FunctionEvaluate(Rmax, array);
            dphi = bspline::Basis::FunctionEvaluate(Rmax, array, 1);

            double k = std::sqrt(2.*energy);


            double phase = phase_shift(k, l, Z, Rmax, phi, dphi);

            phaseshifts[l].push_back({energy, phase});

        }
    }
}


double phase_shift(double k, int l, double Z, double r, maths::complex phi, maths::complex dphi) {
    maths::complex numer = 1.i*phi + dphi / (k + Z/(k*r));      // exp(i*phi)
    maths::complex denom = std::pow(2.0*k*r, 1i*Z/k);

    return std::arg(numer / denom) - k*r + l*maths::Pi/2.0;
 }