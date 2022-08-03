#pragma once

#include "common/objects/potential.h"
#include "common/utility/banded_matrix.h"

#include <functional>
#include <vector>

class YukawaPotential : public tdse::Potential::Register<YukawaPotential> {
    double _Z, _D;
    double _x, _y, _z;
    double _r0, _theta0, _phi0;
    int _expansion_lmax;

    static void BuildExpInvR(int N, double Z, double D, BandedMatrix& invR);
    static void BuildExpInvRR(int N, double Z, double D, BandedMatrix& invRR);
    static void FillBlock( int l1, int m1, int l2, int m2,
                            int lmax, int mmax, int N, int order,
                            const std::vector<int>& Ms,
                            const std::vector<int>& mRows,
                            const BandedMatrix& invR,
                            const BandedMatrix& invRR,
                            std::function<maths::complex(int, int, int, int)> YlmYlm,
                            maths::Matrix m);


public:
    YukawaPotential();

    double operator() (double x, double y, double z) const;
    void FillMatrix(maths::Matrix m, int N, const std::vector<int>& Ms = {0}, const std::vector<int>& mRows = {0});
    void FillMatrixGradX(maths::Matrix m, int N, int lmax, const std::vector<int>& Ms = {0}, const std::vector<int>& mRows = {0});
    void FillMatrixGradY(maths::Matrix m, int N, int lmax, const std::vector<int>& Ms = {0}, const std::vector<int>& mRows = {0});
    void FillMatrixGradZ(maths::Matrix m, int N, int lmax, const std::vector<int>& Ms = {0}, const std::vector<int>& mRows = {0});

    static std::string GetName();
    static tdse::Potential::Ptr_t Create(const nlohmann::json& observable);
    static bool Validate(const nlohmann::json& observable);
};