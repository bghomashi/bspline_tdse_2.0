#include "objects/potentials/coulomb/coulomb_pot.h"
#include "common/utility/index_manip.h"
#include "common/utility/logger.h"
#include "common/utility/spherical_harmonics.h"
#include "common/bspline/bspline.h"

using namespace tdse;
using namespace maths;

void CoulombPotential::FillMatrixGradX(Matrix m, int N, int lmax, const std::vector<int>& Ms, const std::vector<int>& mRows) {
    int order = bspline::Basis::GetOrder();
    // if we get here this is a full 3D calculations
    // so Ms = [-mmax,mmax]
    int mmax = Ms.back();

    if (_isCentral) {
        // the same for each L-M block so cache it
        BandedMatrix invRR(N, 2*order-1);
        BuildInvRR(N, _Z, invRR);

        for (int m1 : Ms) {
            for (int l1 = std::abs(m1); l1 <= lmax; l1++) {
                FillBlock(l1, m1, l1+1, m1+1,
                          lmax, mmax, N, order, 
                          Ms, mRows, invRR, YlmXYlm, m);
                FillBlock(l1, m1, l1-1, m1+1,
                          lmax, mmax, N, order, 
                          Ms, mRows, invRR, YlmXYlm, m);
                FillBlock(l1, m1, l1+1, m1-1,
                          lmax, mmax, N, order, 
                          Ms, mRows, invRR, YlmXYlm, m);
                FillBlock(l1, m1, l1-1, m1-1,
                          lmax, mmax, N, order, 
                          Ms, mRows, invRR, YlmXYlm, m);
            }
            // int m1Block = mRows[i]/N;                     // number of l-blocks to skip (rows)
            // int m2;
            // int m2Block;                                    // number of l-blocks to skip (cols)
            // int blockRow, blockCol, l2;

            // if (m1+1 <= mmax) {                            // then there is a coupling
            //     m2 = m1+1;
            //     m2Block = RowFrom(m2, Ms, mRows)/N;
            //     // for each l-block (block rows)
            //     for (int l1 = std::abs(m1); l1 <= lmax; l1++) {         // m1=m2-1, l1=l2-1
            //         // check one l-block up

            //         if (l1+1 <= lmax && l1+1 >= std::abs(m2)) {
            //             l2 = l1+1;
            //             blockRow = m1Block + (l1-std::abs(m1));
            //             blockCol = m2Block + (l2-std::abs(m2));

            //             m->FillBandedBlock(order-1, N, blockRow, blockCol, 
            //             [=,&invRR](int row, int col) {
            //                 int i = row % N, j = col % N;
            //                 return invRR[i + j*N]*YlmXYlm(l1,m1,l2,m2);
            //             });
            //         }

            //         // check one l-block down
            //         if (l1-1 >= 0 && l1-1 >= std::abs(m2)) {            // m1=m2-1, l1=l2+1
            //             l2 = l1-1;
            //             blockRow = m1Block + (l1-std::abs(m1));
            //             blockCol = m2Block + (l2-std::abs(m2));

            //             m->FillBandedBlock(order-1, N, blockRow, blockCol, 
            //             [=,&invRR](int row, int col) {
            //                 int i = row % N, j = col % N;
            //                 return invRR[i + j*N]*YlmXYlm(l1,m1,l2,m2);
            //             });
            //         }
            //     }
            // }
            // if (m1-1 >= -mmax) {                           // then there is a coupling
            //     m2 = m1-1;
            //     m2Block = RowFrom(m2, Ms, mRows)/N;
            //     // for each l-block (block rows)
            //     for (int l1 = std::abs(m1); l1 <= lmax; l1++) {
            //         // check one l-block up
            //         if (l1+1 <= lmax && l1+1 >= std::abs(m2)) {         // m1=m2+1, l1=l2-1
            //             l2 = l1+1;
            //             blockRow = m1Block + (l1-std::abs(m1));
            //             blockCol = m2Block + (l2-std::abs(m2));

            //             m->FillBandedBlock(order-1, N, blockRow, blockCol, 
            //             [=,&invRR](int row, int col) {
            //                 int i = row % N, j = col % N;
            //                 return invRR[i + j*N]*YlmXYlm(l1,m1,l2,m2);
            //             });
            //         }
            //         // check one l-block down
            //         if (l1-1 >= 0 && l1-1 >= std::abs(m2)) {            // m1=m2+1, l1=l2+1
            //             l2 = l1-1;
            //             blockRow = m1Block + (l1-std::abs(m1));
            //             blockCol = m2Block + (l2-std::abs(m2));

            //             m->FillBandedBlock(order-1, N, blockRow, blockCol, 
            //             [=,&invRR](int row, int col) {
            //                 int i = row % N, j = col % N;
            //                 return invRR[i + j*N]*YlmXYlm(l1,m1,l2,m2);
            //             });
            //         }
            //     }
            
        }

        m->AssembleBegin();
        m->AssembleEnd();
    } else {
        assert(_isCentral && "only supporting central potentials.");
    }
}
