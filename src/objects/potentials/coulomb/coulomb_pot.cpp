#include "objects/potentials/coulomb/coulomb_pot.h"
#include "common/utility/index_manip.h"
#include "common/utility/logger.h"
#include "common/utility/spherical_harmonics.h"
#include "common/bspline/bspline.h"
#include "common/system_state/system_state.h"

#include <iostream>
#include <functional>

using namespace tdse;
using namespace maths;


CoulombPotential::CoulombPotential() : _Z(0), _x(0), _y(0), _z(0), _r0(0), _theta0(0), _phi0(0), _expansion_lmax(0) {}

double CoulombPotential::operator() (double x, double y, double z) const {
    double r = std::sqrt((x-_x)*(x-_x) + (y-_y)*(y-_y) + (z-_z)*(z-_z));
    return -_Z/r;
}

// FIX: maybe this way of thinking about it is not very efficient
void CoulombPotential::FillMatrix(Matrix m, int N, const std::vector<int>& Ms, const std::vector<int>& mRows) {
    int order = bspline::Basis::GetOrder();

    if (_isCentral) {
        // the same for each L-M block so cache it
        BandedMatrix invR(N, 2*order-1);
        BuildInvR(N, _Z, invR);

        // fill across the main diagonal-blocks
        m->FillBandedBlock(order-1, N, [&](int row, int col) {
            int i, j;
            i = row % N;
            j = col % N;
            return invR(i, j);
        });
    } else if (_isAxial) {
        int lmax = SystemState::GetBasisLmax();
        if (_expansion_lmax == -1) _expansion_lmax = lmax;
        std::vector<BandedMatrix> r_matrix(_expansion_lmax+1);
        for (int l = 0; l <= _expansion_lmax; l++) 
            r_matrix[l] = R_matrix(N, order, l, _r0);

        
        for (int m1 : Ms) {                                         // for each m-block (block rows)      
            for (int l1 = std::abs(m1); l1 <= lmax; l1++) {         // for each l1-block (block rows)
                for (int l2 = std::abs(m1); l2 <= lmax; l2++) {     // for each l2-block (block cols)
	                // diagonal in m's - each m-block is different but not coupled
                    int m1Block = RowFrom(m1, Ms, mRows)/N;
                    int blockRow = m1Block + (l1-std::abs(m1));
                    int blockCol = m1Block + (l2-std::abs(m1));

                    m->FillBandedBlock(order-1, N, blockRow, blockCol, [=,&r_matrix](int row, int col) {
                        int i = row % N, j = col % N;
                        return il1m1_jl2m2(i, l1, m1, j, l2, m1, r_matrix);
                    });
                }
            }
        }

        m->AssembleBegin();
        m->AssembleEnd();
    } else {
        assert(_isCentral && _isAxial && "only supporting central and axial potentials.");
    }
}




// utility functions
void CoulombPotential::BuildInvR(int N, double Z, BandedMatrix& invR) {
    int order = bspline::Basis::GetOrder();

    invR = BandedMatrix(N, 2*order-1);
    for (int i = 0; i < N; i++) {
        for (int j = i; j < std::min(N, i+order); j++) {
            invR(i, j) = invR(j, i) = bspline::Basis::Integrate(i+1, j+1, [Z] (complex r) -> complex {
                return -Z/r;
            });
        }
    }
}
void CoulombPotential::BuildInvRR(int N, double Z, BandedMatrix& invRR) {
    int order = bspline::Basis::GetOrder();
    
    invRR = BandedMatrix(N, 2*order-1);
    for (int i = 0; i < N; i++) {
        for (int j = i; j < std::min(N, i+order); j++) {
            invRR(i, j) = invRR(j, i) = bspline::Basis::Integrate(i+1, j+1, [Z] (complex r) -> complex {
                return Z/r/r;
            });
        }
    }
}

void CoulombPotential::FillBlock( int l1, int m1, int l2, int m2,
                int lmax, int mmax, int N, int order,
                const std::vector<int>& Ms,
                const std::vector<int>& mRows,
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
        [=,&invRR](int row, int col) {
            int i = row % N, j = col % N;
            return invRR(i, j)*YlmYlm(l1,m1,l2,m2);
        });
    }
}


std::string CoulombPotential::GetName() {
    return "coulomb";
}
Potential::Ptr_t CoulombPotential::Create(const nlohmann::json& pot_term) {
    auto coulomb_pot = new CoulombPotential();

    coulomb_pot->_expansion_lmax = -1;
    if (pot_term.contains("lmax"))
        coulomb_pot->_expansion_lmax = pot_term["lmax"];

    coulomb_pot->_Z = pot_term["Z"].get<double>();
    if (pot_term.contains("location")) {
        double x = pot_term["location"][0].get<double>();
        double y = pot_term["location"][1].get<double>();
        double z = pot_term["location"][2].get<double>();


        coulomb_pot->_x = x;
        coulomb_pot->_y = y;
        coulomb_pot->_z = z;

        coulomb_pot->_r0 = std::sqrt(x*x + y*y + z*z);
        // check domain of arctan arguements to avoid domain error.
        coulomb_pot->_theta0 = (x == 0 && y == 0 && z == 0 ? 0. : std::atan2(std::sqrt(x*x + y*y), z));
        coulomb_pot->_phi0 = (x == 0 && y == 0 ? 0. : std::atan2(y, x));

        coulomb_pot->_isCentral = coulomb_pot->_x == 0 && coulomb_pot->_y == 0 && coulomb_pot->_z == 0;
        coulomb_pot->_isAxial   = coulomb_pot->_x == 0 && coulomb_pot->_y == 0;
    }

    return Potential::Ptr_t(coulomb_pot);
}
bool CoulombPotential::Validate(const nlohmann::json& pot_term) {
    if (!(pot_term.contains("Z") && pot_term["Z"].is_number())) {
        LOG_CRITICAL("\"coulomb\" potential must contain number entry: Z");
        return false;
    }
    if (pot_term.contains("location") && !pot_term["location"].is_array()) {
        LOG_CRITICAL("Optional entry \"location\" must be a vector in Coloumb term.");
        return false;
    }
    return true;
}
