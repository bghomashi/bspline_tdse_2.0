
#include "objects/propagators/crank_nicolson/crank_nicolson.h"
#include "common/utility/index_manip.h"
#include "common/tdse/simulation.h"
#include "common/bspline/bspline.h"
#include "common/system_state/system_state.h"


using namespace std::complex_literals;
using namespace tdse;
using Maths = maths::Factory;
using namespace maths;

void CrankNicolson::FillU0(Matrix& U0) {
    int N = bspline::Basis::GetNumBSplines();
    int order = bspline::Basis::GetOrder();
    int dof = Simulation::GetDOF();
    int mmax = SystemState::GetBasisMmax();
    int lmax = SystemState::GetBasisLmax();
    auto& Ms = Simulation::GetMs();
    auto& MRows = Simulation::GetMRows();
    auto pol = Simulation::GetPolarization();

    auto fOne = [](int r, int c) -> maths::complex { return 1.; };

    // fill the main diagonal blocks
    U0->FillBandedBlock(order-1, N, fOne);

    if (pol[Z]) {              // if there is z-polarization
        // fill the band above/block the main diagonal
        U0->FillBandedBlock(order-1, N, 1, fOne);
        U0->FillBandedBlock(order-1, N, -1, fOne);;
    }
    if (pol[X] || pol[Y]) {              // if there is x/y-polarization
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

                        U0->FillBandedBlock(order-1, N, blockRow, blockCol, fOne);
                    }

                    // check one l-block down
                    if (l1-1 >= 0 && l1-1 >= std::abs(m2)) {
                        l2 = l1-1;
                        blockRow = m1Block + (l1-std::abs(m1));
                        blockCol = m2Block + (l2-std::abs(m2));

                        U0->FillBandedBlock(order-1, N, blockRow, blockCol, fOne);
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

                        U0->FillBandedBlock(order-1, N, blockRow, blockCol, fOne);
                    }
                    // check one l-block down
                    if (l1-1 >= 0 && l1-1 >= std::abs(m2)) {
                        l2 = l1-1;
                        blockRow = m1Block + (l1-std::abs(m1));
                        blockCol = m2Block + (l2-std::abs(m2));

                        U0->FillBandedBlock(order-1, N, blockRow, blockCol, fOne);
                    }
                }
            }
        }


    }
    U0->AssembleBegin();
    U0->AssembleEnd();

} 