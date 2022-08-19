#include "objects/potentials/yukawa/yukawa_pot.h"
#include "common/utility/index_manip.h"
#include "common/utility/logger.h"
#include "common/utility/spherical_harmonics.h"
#include "common/bspline/bspline.h"

using namespace maths;
using namespace tdse;

YukawaPotential::YukawaPotential() : _Z(0), _D(0), _x(0), _y(0), _z(0), _r0(0), _theta0(0), _phi0(0), _expansion_lmax(10) {}

double YukawaPotential::operator() (double x, double y, double z) const {
    double r = std::sqrt ((x-_x)*(x-_x) + (y-_y)*(y-_y) + (z-_z)*(z-_z));
    return (-_Z/r)*std::exp(-_D*r);                                        // currently only central is supported
}

void YukawaPotential::FillMatrix(Matrix m, int N, const std::vector<int>& Ms, const std::vector<int>& mRows) {
    int order = bspline::Basis::GetOrder();

    if (_isCentral) {
        // the same for each L-M block so cache it
        BandedMatrix invR(N, 2*order-1);
        BuildExpInvR(N, _Z, _D, invR);

        m->FillBandedBlock(order-1, N, [&](int row, int col) {
            int i, j, l1, l2, m1, m2;
            ILMFrom(row, i, l1, m1, N, Ms, mRows);
            ILMFrom(col, j, l2, m2, N, Ms, mRows);
            return invR(i, j);
        });
    } else {
        assert(_isCentral && "only supporting central potentials.");
    }
}



// utility functions
void YukawaPotential::BuildExpInvR(int N, double Z, double D, BandedMatrix& invR) {
    int order = bspline::Basis::GetOrder();

    invR = BandedMatrix(N, 2*order-1);
    for (int i = 0; i < N; i++) {
        for (int j = i; j < std::min(N, i+order); j++) {
            invR(i, j) = invR(j, i) = bspline::Basis::Integrate(i+1, j+1, [Z,D] (complex r) -> complex {
                return (-Z/r)*std::exp(-D*r);
            });
        }
    }
}
void YukawaPotential::BuildExpInvRR(int N, double Z, double D, BandedMatrix& invRR) {
    int order = bspline::Basis::GetOrder();

    invRR = BandedMatrix(N, 2*order-1);
    for (int i = 0; i < N; i++) {
        for (int j = i; j < std::min(N, i+order); j++) {
            invRR(i, j) = invRR(j, i) = bspline::Basis::Integrate(i+1, j+1, [Z,D] (complex r) -> complex {
                return (Z/r/r)*std::exp(-D*r);
            });
        }
    }
}

void YukawaPotential::FillBlock( int l1, int m1, int l2, int m2,
                int lmax, int mmax, int N, int order,
                const std::vector<int>& Ms,
                const std::vector<int>& mRows,
                const BandedMatrix& invR,
                const BandedMatrix& invRR,
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
        [=,&invR,&invRR](int row, int col) {
            int i = row % N, j = col % N;
            return (invR(i, j) + invRR(i, j))*YlmYlm(l1,m1,l2,m2);
        });
    }
}


const std::string YukawaPotential::Name() const {
    return GetName();
}
std::string YukawaPotential::GetName() {
    return "yukawa";
}
Potential::Ptr_t YukawaPotential::Create(const nlohmann::json& pot_term) {
    auto yukawa_pot = new YukawaPotential();

    yukawa_pot->_D = pot_term["decay"].get<double>();
    yukawa_pot->_Z = pot_term["Z"].get<double>();
    if (pot_term.contains("location")) {
        // double x = pot_term["location"][0].get<double>();
        // double y = pot_term["location"][1].get<double>();
        // double z = pot_term["location"][2].get<double>();

        // yukawa_pot->_x = x;
        // yukawa_pot->_y = y;
        // yukawa_pot->_z = z;

        // yukawa_pot->_r0 = std::sqrt(x*x + y*y + z*z);
        // // check domain of arctan arguements to avoid domain error.
        // yukawa_pot->_theta0 = (x == 0 && y == 0 && z == 0 ? 0. : std::atan2(std::sqrt(x*x + y*y), z));
        // yukawa_pot->_phi0 = (x == 0 && y == 0 ? 0. : std::atan2(y, x));

        // yukawa_pot->_isCentral = yukawa_pot->_x == 0 && yukawa_pot->_y == 0 && yukawa_pot->_z && 0;
        // yukawa_pot->_isAxial   = yukawa_pot->_x == 0 && yukawa_pot->_y == 0;
    }

    return Potential::Ptr_t(yukawa_pot);
}
bool YukawaPotential::Validate(const nlohmann::json& pot_term) {
    if (!(pot_term.contains("Z") && pot_term["Z"].is_number())) {
        LOG_CRITICAL("\"yukawa\" observable must contain number entry: Z");
        return false;
    }
    if (!pot_term.contains("decay") || !pot_term["decay"].is_number() || pot_term["decay"].get<int>() < 0) {
        LOG_CRITICAL("\"yukawa\" potential must contain (positive) number entry: decay");
        return false;
    }
    // if (pot_term.contains("location") && !pot_term["location"].is_array()) {
    //     LOG_CRITICAL("Optional entry \"location\" must be a 3-vector in \"yukawa\" term.");
    //     return false;
    // }
    return true;
}
