#include "objects/potentials/coulomb/coulomb_pot.h"
#include "common/utility/index_manip.h"
#include "common/utility/logger.h"
#include "common/utility/spherical_harmonics.h"
#include "common/bspline/bspline.h"

using namespace tdse;
using namespace maths;


void CoulombPotential::FillMatrixGradY(Matrix m, int N, int lmax, const std::vector<int>& Ms, const std::vector<int>& mRows) {    
    int order = bspline::Basis::GetOrder();
    int mmax = Ms.back();

    if (_isCentral) {
        // the same for each L-M block so cache it
        BandedMatrix invRR(N, 2*order-1);
        BuildInvRR(N, _Z, invRR);

        // for each m-block (block rows)
        for (int m1 : Ms) {
            for (int l1 = std::abs(m1); l1 <= lmax; l1++) {
                FillBlock(l1, m1, l1+1, m1+1,
                          lmax, mmax, N, order, 
                          Ms, mRows, invRR, YlmYYlm, m);
                FillBlock(l1, m1, l1-1, m1+1,
                          lmax, mmax, N, order, 
                          Ms, mRows, invRR, YlmYYlm, m);
                FillBlock(l1, m1, l1+1, m1-1,
                          lmax, mmax, N, order, 
                          Ms, mRows, invRR, YlmYYlm, m);
                FillBlock(l1, m1, l1-1, m1-1,
                          lmax, mmax, N, order, 
                          Ms, mRows, invRR, YlmYYlm, m);
            }
        }

        m->AssembleBegin();
        m->AssembleEnd();
    } else {
        assert(_isCentral && "only supporting central potentials.");
    }
}
