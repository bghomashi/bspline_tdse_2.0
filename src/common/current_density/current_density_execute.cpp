#include "common/current_density/current_density.h"
#include "common/utility/logger.h"
#include "common/utility/index_manip.h"
#include "common/utility/banded_matrix.h"
#include "common/utility/spherical_harmonics.h"
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

// void BuildOverlap(Matrix S);

void CurrentDensity::Execute() {
    int potSymmetry = Simulation::GetPotentialSymmetry();
    if (potSymmetry == Symmetry::Central)
        PopulationsCentral();
    // else if (potSymmetry == Symmetry::Axial)
    //     PopulationsAxial();
}

void CurrentDensity::PopulationsCentral() {
    complex pop;
    std::stringstream ss;
    io::Binary inFile, outFile;
    // ---------- grab a bunch of state variables ------

    // ------------- tabulate m-block rows ---------------
    // count degrees of freedom and 
    // generate row table for m's
    


    // --------- build overlap -------------
    // BuildOverlap(S);

    // ----------- open text file
    int it = 10800;
    if ((inFile = io::Factory::OpenBinary("density_"+std::to_string(it)+".bin", 'r')) == nullptr) {
        LOG_INFO("Failed to open file.");
        return;
    }
    if ((outFile = io::Factory::OpenBinary("current_density_"+std::to_string(it)+".bin", 'w')) == nullptr) {
        LOG_INFO("Failed to open file.");
        return;
    }

    int numGrid;
    double dx, xmin, xmax;
    std::vector<complex> amplitudes;
    std::vector<double> current_x;
    std::vector<double> current_y;

    inFile->Read(&numGrid, sizeof(int)); 
    inFile->Read(&dx, sizeof(double)); 
    inFile->Read(&xmin, sizeof(double)); 
    inFile->Read(&xmax, sizeof(double)); 
    amplitudes.resize(numGrid*numGrid);
    current_x.resize(numGrid*numGrid);
    current_y.resize(numGrid*numGrid);

    for (int i = 0; i < numGrid; i++) {
        for (int j = 0; j < numGrid; j++) {
            // for (int k = 0; k < numGrid; k++)
            {
                double x = xmin + i*dx;
                double y = xmin + j*dx;
                double z = 0;

                inFile->Read(&amplitudes[j+i*numGrid], sizeof(complex)); 
            }            
        }  
    }

    outFile->Write(&numGrid, sizeof(int)); 
    outFile->Write(&dx, sizeof(double)); 
    outFile->Write(&xmin, sizeof(double)); 
    outFile->Write(&xmax, sizeof(double)); 
    complex ddx, ddy;
    double djx, djy;
    for (int i = 0; i < numGrid; i++) {
        for (int j = 0; j < numGrid; j++) {
            // differentiate along x
            if (i == 0)
                ddx = (amplitudes[j+1*numGrid] - amplitudes[j+0*numGrid]) / dx;
            else if (i == numGrid-1)
                ddx = (amplitudes[j+(numGrid-1)*numGrid] - amplitudes[j+(numGrid-2)*numGrid]) / dx;
            else
                ddx = (amplitudes[j+(i+1)*numGrid] - amplitudes[j+(i-1)*numGrid]) / (2*dx);

            // differentiate along y
            if (j == 0)
                ddy = (amplitudes[1+i*numGrid] - amplitudes[0+i*numGrid]) / dx;
            else if (j == numGrid-1)
                ddy = (amplitudes[(numGrid-1)+i*numGrid] - amplitudes[(numGrid-2)+i*numGrid]) / dx;
            else
                ddy = (amplitudes[(j+1)+i*numGrid] - amplitudes[(j-1)+i*numGrid]) / (2*dx);

            djx = std::imag(ddx*amplitudes[j+i*numGrid]);
            djy = std::imag(ddy*amplitudes[j+i*numGrid]);
            for (auto& p : _pulses) {
                djx -= p->A(it*_time_step).x*std::abs(amplitudes[j+i*numGrid])*std::abs(amplitudes[j+i*numGrid]);
                djy -= p->A(it*_time_step).y*std::abs(amplitudes[j+i*numGrid])*std::abs(amplitudes[j+i*numGrid]);
            }
            outFile->Write(&djx, sizeof(double)); 
            outFile->Write(&djy, sizeof(double)); 
        }
    }
    

    
    
    inFile = nullptr;
    outFile = nullptr;
}


// void BuildOverlap(Matrix S) {
//     int N = bspline::Basis::GetNumBSplines();
//     int order = bspline::Basis::GetNumBSplines();
//     BandedMatrix overlapStore(N, 2*order-1);

//     for (int i = 0; i < N; i++) {
//         for (int j = i; j < std::min(N, i+order); j++) {
//             overlapStore(j,i) = overlapStore(i,j) = bspline::Basis::Integrate(i+1, j+1);
//         }
//     }
//     S->FillBandedBlock(order-1, N, [=,&overlapStore](int row, int col) {
//         int i, j;
//         i = row % N;  
//         j = col % N;  

//         return overlapStore(i, j);
//     });
// }

    // BandedMatrix S(N, 2*order-1);
    // BandedMatrix dr(N, 2*order-1);
    // BandedMatrix r1(N, 2*order-1);

    // for (int i = 0; i < N; i++) {
    //     for (int j = i; j < std::min(N, i+order); j++) {
    //         S(j,i) = S(i,j) = bspline::Basis::Integrate(i+1, j+1);
    //         dr(i,j) = bspline::Basis::Integrate(i+1, j+1, 0, 1);
    //         dr(j,i) = bspline::Basis::Integrate(j+1, i+1, 0, 1);
    //         r1(i,j) = r1(j,i) = bspline::Basis::Integrate(j+1, i+1, 0, 0, [](complex x) {
    //             return 1./x;
    //         });
    //     }
    // }


    // for each x {
    //     for each y {
    //         double r = x*x + y*y;
    //         double t = 0;
    //         double p = atan2(y,x);
            

    //         for (int i = 0; i < nodes; i++) {
    //             for (int m1 = -mmin; m1 < mmax; m1++) {
    //                 for (int l1 = std::abs(m1); l1 <= lmax; l1++) {


    //                     for (int j = 0; j < nodes; j++) {
    //                         for (int m2 = -mmin; m2 < mmax; m2++) {
    //                             for (int l2 = std::abs(m2); l2 <= lmax; l2++) {
    //                                 c(i,l1,m1)*c(j,l3,m3)*
    //                                 std::conj(B(i,r))*B(k,r) 
    //                                 Ylm(l1,m1,t,p)*std::conj(Ylm(l2,m2,t,p))*
    //                                 (dr(k,j) + 0.5*(l3*(l3+1) - l2*(l2+1))*r1(k,j))*
    //                                 YlmzYlm(l2,m2,l3,m3);
    //                             }
    //                         }
    //                     }
    //                 }
    //             }
    //         }
    //     }
    // }


