
#include "objects/propagators/crank_nicolson/crank_nicolson.h"
#include "common/utility/index_manip.h"
#include "common/utility/spherical_harmonics.h"
#include "common/utility/logger.h"
#include "common/utility/banded_matrix.h"
#include "common/tdse/simulation.h"
#include "common/bspline/bspline.h"
#include "common/system_state/system_state.h"

using namespace std::complex_literals;
using namespace tdse;
using Maths = maths::Factory;
using namespace maths;

void CrankNicolson::FillInteractionY(Matrix& HI) {
    int N = bspline::Basis::GetNumBSplines();
    int order = bspline::Basis::GetOrder();
    int dof = Simulation::GetDOF();
    int mmax = SystemState::GetBasisMmax();
    int lmax = SystemState::GetBasisLmax();
    auto& Ms = Simulation::GetMs();
    auto& MRows = Simulation::GetMRows();
    auto& potentials = Simulation::GetPotentials();

    // if we made it here we assume there IS m->m+1 coupling
    // so a full -mmax to mmax matrix
    BandedMatrix ddr(N, 2*order-1);
    BandedMatrix invR(N, 2*order-1);

    for (int r = 0; r < N; r++) {
        for (int c = 0; c < N; c++) {
            ddr(r, c) = bspline::Basis::Integrate(r+1, c+1, 0, 1);
            ddr(c, r) = bspline::Basis::Integrate(c+1, r+1, 0, 1);
            invR(r, c) = invR(c, r) = bspline::Basis::Integrate(r+1, c+1, [](complex x) {
                return 1./x;
            });
        }
	}

    // for each m-block (block rows)
    for (int i = 0; i < Ms.size(); i++) {              
        int m1 = Ms[i];
        int m1Block = MRows[i]/N;                     // number of l-blocks to skip (rows)
        int m2;
        int m2Block;                                    // number of l-blocks to skip (cols)
        int blockRow, blockCol, l2;

        if (m1+1 <= mmax) {                            // then there is a coupling
            m2 = m1+1;
            m2Block = RowFrom(m2, Ms, MRows)/N;
            // for each l-block (block rows)
            for (int l1 = std::abs(m1); l1 <= lmax; l1++) {
                // check one l-block up
                if (l1+1 <= lmax && l1+1 >= std::abs(m2)) {
                    l2 = l1+1;
                    blockRow = m1Block + (l1-std::abs(m1));
                    blockCol = m2Block + (l2-std::abs(m2));

                    HI->FillBandedBlock(order-1, N, blockRow, blockCol, 
                    [=,&ddr,&invR](int row, int col) {
                        int i = row % N, j = col % N;
                        
                        // m1=m2-1, l1=l2-1
                        double a = -sqrt((l2+m2) * (l2+m2-1) / (2.*l2 + 1.) / (2.*l2 - 1.));

                        return (ddr(i, j) + double(l2)*invR(i, j))*YlmYYlm(l1,m1,l2,m2);
                    });
                }

                // check one l-block down
                if (l1-1 >= 0 && l1-1 >= std::abs(m2)) {
                    l2 = l1-1;
                    blockRow = m1Block + (l1-std::abs(m1));
                    blockCol = m2Block + (l2-std::abs(m2));

                    HI->FillBandedBlock(order-1, N, blockRow, blockCol, 
                    [=,&ddr,&invR](int row, int col) {
                        int i = row % N, j = col % N;
                        
                        // m1=m2-1, l1=l2+1
                        double a = sqrt((l2-m2+1)*(l2-m2+2) / (2.*l2 + 1.) / (2.*l2 + 3.));

                        return (ddr(i, j) - double(l2+1.)*invR(i, j))*YlmYYlm(l1,m1,l2,m2);
                    });
                }
            }
        }
        if (m1-1 >= -mmax) {                           // then there is a coupling
            m2 = m1-1;
            m2Block = RowFrom(m2, Ms, MRows)/N;
            // for each l-block (block rows)
            for (int l1 = std::abs(m1); l1 <= lmax; l1++) {
                // check one l-block up
                if (l1+1 <= lmax && l1+1 >= std::abs(m2)) {
                    l2 = l1+1;
                    blockRow = m1Block + (l1-std::abs(m1));
                    blockCol = m2Block + (l2-std::abs(m2));

                    HI->FillBandedBlock(order-1, N, blockRow, blockCol, 
                    [=,&ddr,&invR](int row, int col) {
                        int i = row % N, j = col % N;
                        
                        // m1=m2+1, l1=l2-1
                        double a = sqrt((l2-m2)*(l2-m2-1) / (2.*l2 + 1.) / (2.*l2 - 1.));

                        return (ddr(i, j) + double(l2)*invR(i, j))*YlmYYlm(l1,m1,l2,m2);
                    });
                }
                // check one l-block down
                if (l1-1 >= 0 && l1-1 >= std::abs(m2)) {
                    l2 = l1-1;
                    blockRow = m1Block + (l1-std::abs(m1));
                    blockCol = m2Block + (l2-std::abs(m2));

                    HI->FillBandedBlock(order-1, N, blockRow, blockCol, 
                    [=,&ddr,&invR](int row, int col) {
                        int i = row % N, j = col % N;
                        
                        // m1=m2+1, l1=l2+1
                        double a = -sqrt((l2+m2+1)*(l2+m2+2) / (2.*l2 + 1.) / (2.*l2 + 3.));

                        return (ddr(i, j) - double(l2+1.)*invR(i, j))*YlmYYlm(l1,m1,l2,m2);
                    });
                }
            }
        }
    }

    HI->AssembleBegin();
    HI->AssembleEnd();
}
