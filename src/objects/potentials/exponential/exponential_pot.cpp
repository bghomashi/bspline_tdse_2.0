#include "objects/potentials/exponential/exponential_pot.h"
#include "common/utility/index_manip.h"
#include "common/utility/logger.h"
#include "common/utility/spherical_harmonics.h"
#include "common/bspline/bspline.h"

using namespace tdse;
using namespace maths;


ExponentialPotential::ExponentialPotential() : _Z(0), _D(0)/*, _x(0), _y(0), _z(0)*/ {}


double ExponentialPotential::operator() (double x, double y, double z) const {
    double r = std::sqrt (x*x + y+y + z*z);
    return -_Z*std::exp(-_D*r);
}
void ExponentialPotential::FillMatrix(Matrix m, int N, const std::vector<int>& Ms, const std::vector<int>& mRows) {
    int order = bspline::Basis::GetOrder();

    if (_isCentral) {
        // the same for each L-M block so cache it
        BandedMatrix expR(N, 2*order-1);
        BuildExpR(N, _Z, _D, expR);

        m->FillBandedBlock(order-1, N, [&](int row, int col) {
            int i, j, l1, l2, m1, m2;
            ILMFrom(row, i, l1, m1, N, Ms, mRows);
            ILMFrom(col, j, l2, m2, N, Ms, mRows);
            return expR(i, j);
        });
    } else {
        assert(_isCentral && "only supporting central potentials.");
    }
}



// utility functions
void ExponentialPotential::BuildExpR(int N, double Z, double D, BandedMatrix& expR) {
    int order = bspline::Basis::GetOrder();
    
    for (int i = 0; i < N; i++) {
        for (int j = i; j < std::min(N, i+order); j++) {
            expR(i, j) = expR(j, i) = bspline::Basis::Integrate(i+1, j+1, [Z,D] (complex r) -> complex {
                return -Z*std::exp(-D*r);
            });
        }
    }
}

void ExponentialPotential::FillBlock( int l1, int m1, int l2, int m2,
                int lmax, int mmax, int N, int order,
                const std::vector<int>& Ms,
                const std::vector<int>& mRows,
                const BandedMatrix& expR,
                std::function<complex(int, int, int, int)> YlmYlm,
                Matrix m) {
    if (l2 <= lmax && l2 >= std::abs(m2) &&
        l2 >= 0 && l2 >= std::abs(m2) &&
        m2 <= mmax && m2 >= -mmax) {
        int m1Block = RowFrom(m1, Ms, mRows)/N;
        int m2Block = RowFrom(m2, Ms, mRows)/N;
        int blockRow = m1Block + (l1-std::abs(m1));
        int blockCol = m2Block + (l2-std::abs(m2));

        m->FillBandedBlock(order-1, N, blockRow, blockCol, 
        [=,&expR](int row, int col) {
            int i = row % N, j = col % N;
            return expR(i, j)*YlmYlm(l1,m1,l2,m2);
        });
    }
}


const std::string ExponentialPotential::Name() const {
    return GetName();
}
std::string ExponentialPotential::GetName() {
    return "exponential";
}
Potential::Ptr_t ExponentialPotential::Create(const nlohmann::json& pot_term) {
    auto exp_pot = new ExponentialPotential();
    
    exp_pot->_Z = pot_term["amplitude"].get<double>();
    exp_pot->_D = pot_term["decay"].get<double>();

    return Potential::Ptr_t(exp_pot);
}
bool ExponentialPotential::Validate(const nlohmann::json& pot_term) {
    if (!(pot_term.contains("amplitude") && pot_term["amplitude"].is_number())) {
        LOG_CRITICAL("\"exponential\" potential must contain number entry: amplitude");
        return false;
    }
    if (!pot_term.contains("decay") || !pot_term["decay"].is_number() || pot_term["decay"].get<int>() < 0) {
        LOG_CRITICAL("\"exponential\" potential must contain (positive) number entry: decay");
        return false;
    }
    return true;
}
