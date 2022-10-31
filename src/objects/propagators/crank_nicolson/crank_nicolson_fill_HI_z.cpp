
#include "objects/propagators/crank_nicolson/crank_nicolson.h"
#include "common/utility/index_manip.h"
#include "common/utility/logger.h"
#include "common/tdse/simulation.h"
#include "common/bspline/bspline.h"
#include "common/system_state/system_state.h"
#include "common/utility/banded_matrix.h"

#include "petsc/maths/petsc_common.h"

using namespace std::complex_literals;
using namespace tdse;
using Maths = maths::Factory;
using namespace maths;


void CrankNicolson::FillInteractionZ(Matrix& HI) {
    int N = bspline::Basis::GetNumBSplines();
    int order = bspline::Basis::GetOrder();
    int mmax = SystemState::GetBasisMmax();
    int lmax = SystemState::GetBasisLmax();
    int dof = Simulation::GetDOF();
    auto& Ms = Simulation::GetMs();
    auto& MRows = Simulation::GetMRows();

    BandedMatrix ddr(N, 2*order-1), invR(N, 2*order-1);
    
    // ----------------- cache some common matrix elements -----------
    for (int r = 0; r < N; r++) {
        for (int c = r; c < std::min(N, r+order); c++) {
            ddr(r, c)  = bspline::Basis::Integrate(r+1, c+1, 0, 1);           // <Bi|d/dr|Bj>
            ddr(c, r)  = bspline::Basis::Integrate(c+1, r+1, 0, 1);           // <Bi|d/dr|Bj>
            invR(r, c) = invR(c, r) = bspline::Basis::Integrate(r+1, c+1, [](complex x) {  // <Bi|1/r|Bj>
                return 1./x;
            });
        }
	}

    for (int i = 0; i < Ms.size(); i++) {              
        int m1 = Ms[i];
        int m1Block = MRows[i]/N;                     // number of l-blocks to skip (rows)
        int m2 = m1;
        int m2Block = m1Block;                          // number of l-blocks to skip (cols)
        int blockRow, blockCol, l2;

        for (int l1 = std::abs(m1); l1 <= lmax; l1++) {
            // check one l-block up
            if (l1+1 <= lmax && l1+1 >= std::abs(m2)) {
                l2 = l1+1;
                blockRow = m1Block + (l1-std::abs(m1));
                blockCol = m2Block + (l2-std::abs(m2));

                HI->FillBandedBlock(order-1, N, blockRow, blockCol, 
                [=,&ddr,&invR](int row, int col) {
                    int i = row % N, j = col % N;
                    
                    // m1=m2, l1=l2-1
                    double a = sqrt((l2+m2) * (l2-m2) / (2.*l2 + 1.) / (2.*l2 - 1.));

                    return (ddr(i, j) + double(l2)*invR(i, j))*a;
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
                    // m1=m2, l1=l2+1
                    double a = sqrt((l2+m2+1.)*(l2 - m2 + 1.) / (2.*l2 + 1.) / (2.*l2 + 3.));

                    return (ddr(i, j) - double(l2+1)*invR(i, j))*a;
                });
            }
        }
    }

    HI->AssembleBegin();
    HI->AssembleEnd();
    
    // MatView(std::dynamic_pointer_cast<PetscMatrix>(HI)->_petsc_mat, 0); exit(0);
}
