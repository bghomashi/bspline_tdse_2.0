#ifndef __COULOMB_H__
#define __COULOMB_H__

#include "common/objects/potential.h"
#include "common/utility/banded_matrix.h"
#include <functional>
#include <vector>

class CoulombPotential : public tdse::Potential::Register<CoulombPotential> {
    double _Z;
    double _x, _y, _z;
    double _r0, _theta0, _phi0;
    int _expansion_lmax;

    static void BuildInvR(int N, double Z, BandedMatrix& invR);
    static void BuildInvRR(int N, double Z, BandedMatrix& invRR);
    static void FillBlock(   
                            int l1, int m1, int l2, int m2,
                            int lmax, int mmax, int N, int order,
                            const std::vector<int>& Ms,
                            const std::vector<int>& mRows,
                            const BandedMatrix& invRR,
                            std::function<maths::complex(int, int, int, int)> YlmYlm,
                            maths::Matrix m);
                    
    maths::complex R_element(int i, int j, int l, double R0);
    maths::complex il1m1_jl2m2(int i, int l1, int m1, int j, int l2, int m2, const std::vector<BandedMatrix>& rmatrix);
    BandedMatrix R_matrix(int N, int order, int l, double R0);
    maths::Matrix PotentialMatrixCylindrical(int m1);
public:
    CoulombPotential();

    inline double Z() const {
        return _Z;
    }
    double operator() (double x, double y, double z) const;
    void FillMatrix(maths::Matrix m, int N, const std::vector<int>& Ms = {0}, const std::vector<int>& mRows = {0});
    void FillMatrixGradX(maths::Matrix m, int N, int lmax, const std::vector<int>& Ms = {0}, const std::vector<int>& mRows = {0});
    void FillMatrixGradY(maths::Matrix m, int N, int lmax, const std::vector<int>& Ms = {0}, const std::vector<int>& mRows = {0});
    void FillMatrixGradZ(maths::Matrix m, int N, int lmax, const std::vector<int>& Ms = {0}, const std::vector<int>& mRows = {0});

    const std::string Name() const;
    static std::string GetName();
    static tdse::Potential::Ptr_t Create(const nlohmann::json& observable);
    static bool Validate(const nlohmann::json& observable);
};

#endif