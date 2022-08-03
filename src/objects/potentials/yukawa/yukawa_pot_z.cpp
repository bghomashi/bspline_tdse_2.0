#include "objects/potentials/yukawa/yukawa_pot.h"
#include "common/utility/index_manip.h"
#include "common/utility/logger.h"
#include "common/utility/spherical_harmonics.h"
#include "common/bspline/bspline.h"

using namespace maths;
using namespace tdse;


void YukawaPotential::FillMatrixGradZ(Matrix m, int N, int lmax, const std::vector<int>& Ms, const std::vector<int>& mRows) {
    int order = bspline::Basis::GetOrder();
    int mmax = Ms.back();

    if (_isCentral) {
        // the same for each L-M block so cache it
        BandedMatrix invR, invRR;
        BuildExpInvR(N, _Z, _D, invR);
        invR *=-_D;                 // from product rule
        BuildExpInvRR(N, _Z, _D, invRR);


        for (int m1 : Ms) {
            for (int l1 = std::abs(m1); l1 <= lmax; l1++) { 
                FillBlock(l1, m1, l1+1, m1,
                          lmax, mmax, N, order, 
                          Ms, mRows, invR, invRR, 
                          YlmZYlm, m);
                FillBlock(l1, m1, l1-1, m1,
                          lmax, mmax, N, order, 
                          Ms, mRows, invR, invRR, 
                          YlmZYlm, m);
            }
        }
        m->AssembleBegin();
        m->AssembleEnd();
    } else {
        assert(_isCentral && "only supporting central potentials.");
    } 
}
