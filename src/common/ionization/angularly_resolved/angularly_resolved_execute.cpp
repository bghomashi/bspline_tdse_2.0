#include "common/ionization/ionization.h"
#include "common/ionization/angularly_resolved/angularly_resolved.h"
#include "common/maths/math_common.h"
#include "common/file_io/io_factory.h"
#include "common/utility/spherical_harmonics.h"
#include "common/bspline/bspline.h"
#include "common/utility/logger.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>

using namespace std::complex_literals;

// ---------- declare some utility functions ---------------
template <typename T>
T Lerp(double x, double x0, double x1, T y0, T y1) {
    T m = (y1 - y0) / (x1 - x0);
    return y0 + m*(x - x0);
}
int FindInterval(double E, const std::vector<Population>& vals);
int FindInterval(double E, const std::vector<std::pair<double,double>>& vals);
maths::complex InterpolateAmplitude(double E, const Population& a, const Population& b);
double InterpolatePhase(double E, const std::pair<double, double>& a, const std::pair<double, double>& b);
maths::complex pe_amp(  double k, double theta, double phi, 
                        int lmax, int mmax, 
                        std::unordered_map<int, std::vector<std::pair<double, double>>>& d, 
                        std::unordered_map<int, std::unordered_map<int, std::vector<Population>>>& amp);

//---------------------------------------------------------


void AngularlyResolved::Execute() {
    auto& populations = _parent->populations;
    auto& phaseshifts = _parent->phaseshifts;
    int lmax = _parent->lmax;
    int mmax = _parent->mmax;
    int nx, nz;
    std::vector<double>  xs;
    std::vector<double>  zs;
    std::vector<std::vector<double>> ionization_out;

    // --------------- x grid ------------
    nx = (_xmax - _xmin) / _xstep + 1;
    xs.resize(nx);
    for (int i = 0; i < xs.size(); i++)
        xs[i] = _xmin + _xstep*i;

    // --------------- y grid ------------
    nz = (_zmax - _zmin) / _zstep + 1;
    zs.resize(nz);
    for (int i = 0; i < zs.size(); i++)
        zs[i] = _zmin + _zstep*i;

    // -------------- output array --------
    ionization_out.resize(xs.size());
    for (auto& io : ionization_out)
        io.resize(zs.size(), 0);


    for (int j = 0; j < zs.size(); j++) {
        for (int i = 0; i < xs.size(); i++) {
            double px = xs[i];
            double pz = zs[j];

            double phi = (px != 0 ? 
                         // atan2(0, px) : 0);
                         atan2(px, 0) : 0);
            double theta = (px != 0 || pz != 0 ? 
                            atan2(std::abs(px), pz) : 0);
            double k = sqrt(px*px + pz*pz);

            ionization_out[i][j] = std::abs(pe_amp(k, theta, phi, 
                lmax, mmax, 
                phaseshifts, 
                populations));

            ionization_out[i][j] *= ionization_out[i][j];
        }
    }


    std::ofstream outfile(_output_filename);
    outfile << std::setprecision(8) << std::scientific;
    
    for (int i = 0; i < xs.size(); i++) {
        double px = xs[i];
        for (int j = 0; j < zs.size(); j++) {
            double pz = zs[j];
            double pop = ionization_out[i][j];

            outfile << px << "\t" << pz << "\t" << pop << std::endl;
        }
        outfile << std::endl;
    }

}

// ----------- define utility functions -------------------



maths::complex pe_amp_lm(int l, int m, double k, 
    double theta, double phi, 
    double dkl, maths::complex amp_klm) {
        if (k == 0) return 0;
    // maths::complex phase =  std::exp(1.i*dkl);/*std::exp(0.5i*(maths::Pi*l))*/;    // spencer thinks this needs to be conjugated
    maths::complex phase =  std::exp(0.5i*(maths::Pi*l));
    maths::complex result = Ylm(l, m, theta, phi)*amp_klm*phase/k;
    // maths::complex result = dkl;
    // maths::complex result = amp_klm;
    return result;
}

maths::complex pe_amp(  double k, double theta, double phi, 
                        int lmax, int mmax, 
                        std::unordered_map<int, std::vector<std::pair<double, double>>>& d, 
                        std::unordered_map<int, 
                        std::unordered_map<int, 
                        std::vector<Population>>>& amp) {
    maths::complex total_amp(0);


    for (int m = -mmax; m <= mmax; m++) {
        for (int l = std::abs(m); l <= lmax; l++) {
            int amp_interval = FindInterval(0.5*k*k, amp[l][m]);
            int phase_interval = FindInterval(0.5*k*k, d[l]);

            // if (amp_interval == -1)
            //     continue;
                

            // defaults
            maths::complex amp_klm = amp[l][m].back().amplitude;
            double phase = d[l].back().second;

            if (amp_interval == -1) 
                amp_klm = InterpolateAmplitude(0.5*k*k, Population{0,0,0}, amp[l][m][0]);
            else if (amp_interval < amp[l][m].size())
                amp_klm = InterpolateAmplitude(0.5*k*k, amp[l][m][amp_interval], amp[l][m][amp_interval+1]);
            else 
                LOG_CRITICAL("not enough eigenstates. Energy : " + std::to_string(0.5*k*k));
            
            if (amp_interval == -1)
                phase = InterpolatePhase(0.5*k*k, 
                {0,0}, 
                {amp[l][m][amp_interval].energy, amp[l][m][amp_interval+1].phaseshift});
            else if (amp_interval < d[l].size())
                phase = InterpolatePhase(0.5*k*k, 
                {amp[l][m][amp_interval].energy, amp[l][m][amp_interval].phaseshift}, 
                {amp[l][m][amp_interval+1].energy, amp[l][m][amp_interval+1].phaseshift});
                // phase = InterpolatePhase(0.5*k*k, d[l][phase_interval], d[l][phase_interval+1]);

            total_amp += pe_amp_lm(l, m, k, theta, phi, phase, amp_klm);
        }
    }

    return total_amp;
}


maths::complex InterpolateAmplitude(double E, const Population& a, const Population& b) {
    /*
    LOG_CRITICAL(std::string("Lerp amplitude. : (") 
        + std::to_string(E) + ", " + std::to_string(
        std::abs(Lerp(E, a.energy, b.energy, a.amplitude, b.amplitude))) + ") (" + 
        std::to_string(a.energy) + ", " + std::to_string(std::abs(a.amplitude)) + ") (" + 
        std::to_string(b.energy) + ", " + std::to_string(std::abs(b.amplitude)) + ")"
    );
    */
    auto m = (b.amplitude - a.amplitude) / (b.energy - a.energy); //DX + i*DY
    auto result = a.amplitude + m*(E -  a.energy);
    return result;
    // return Lerp(E, a.energy, b.energy, a.amplitude, b.amplitude);
}

double InterpolatePhase(double E, const std::pair<double, double>& a, const std::pair<double, double>& b) {
    return Lerp(E, a.first, b.first, a.second, b.second);
}

int FindInterval(double E, const std::vector<Population>& vals) {
    for (int i = 0; i < vals.size(); i++) {
        if (vals[i].energy > E)
            return i-1;
    }
    return vals.size();
}
int FindInterval(double E, const std::vector<std::pair<double,double>>& vals) {
    for (int i = 0; i < vals.size(); i++) {
        if (vals[i].first > E)
            return i-1;
    }
    return vals.size();
}
