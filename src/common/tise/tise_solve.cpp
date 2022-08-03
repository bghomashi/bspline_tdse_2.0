
#include "common/tise/tise.h"
#include "common/utility/banded_matrix.h"
#include "common/utility/logger.h"
#include "common/utility/profiler.h"
#include "common/utility/file_exists.h"
#include "common/tdse/simulation.h"
#include "common/file_io/io_factory.h"
#include "common/maths/math_algebra.h"
#include "common/maths/eigen_solver.h"
#include "common/bspline/bspline.h"
#include "common/system_state/system_state.h"

#include <string>
#include <sstream>
#include <complex>
#include <iomanip>

using namespace tdse;
using namespace maths;

void tise::TISE::Solve() {
    int N = bspline::Basis::GetNumBSplines();
    int order = bspline::Basis::GetOrder();
    int lmax = SystemState::GetBasisLmax();
    int numNstates = SystemState::GetEigenStateBoundNmax();
    int numLstates = SystemState::GetEigenStateBoundLmax();

    Profile::Push("TISE::Solve");

    // get some info about the most general potential in the list
    bool isCentral = true, isAxial = true;
    for (auto& p : _potentials) {
        if (!p->isCentral())
            isCentral = false;
        else if (!p->isAxial())
            isAxial = false;
    }

    // ---------------- estimating memory ---------------
    // only looking at the matrices since the small variables
    // are insignificant.
    int maxNumOfEigenState = (2*numNstates - numLstates)*(numLstates + 1)/2; // sum_{l=0}^Lmax (Nmax-l)
    int bandsPerBlock = 2*order-1;
    int numBlocks = isCentral ? 1 : lmax + 1;
    int totalCols = bandsPerBlock*numBlocks;
    int totalRows = N*numBlocks;
    int totalNumPerMatrix = totalRows*totalCols;
    int totalFullMatrices = 2 + (isCentral ? _potentials.size() : 1);
    int totalBlockCaches = (isCentral ? 2 : 3);

    double memory = totalFullMatrices*totalNumPerMatrix   + // full matrices
                    totalBlockCaches*N*bandsPerBlock +      // cached blocks
                    totalRows*maxNumOfEigenState +          // eigenvectors
                    maxNumOfEigenState;                     // eigenvalues
    memory = memory*sizeof(complex)/1024./1024./1024.;      // gigabytes
    
    LOG_INFO("Estimated memory required: " + std::to_string(memory) + " GB.");

    if (isCentral)
        SolveCentral();
    else if (isAxial)
        SolveAxial();
    else
        assert((isCentral || isAxial) && "Only central and axial potentials are supported.");
    
    Profile::Pop("TISE::Solve");
}
