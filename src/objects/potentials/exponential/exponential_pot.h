#pragma once

#include "common/objects/potential.h"
#include "common/maths/math_common.h"
#include "common/utility/banded_matrix.h"

#include <functional>
#include <vector>

// currently this is centered at the origin
// - probably add arbitrary shift at some point
class ExponentialPotential : public tdse::Potential::Register<ExponentialPotential> {
    double _Z, _D;

    static void BuildExpR(int N, double Z, double D, BandedMatrix& expR);
    static void FillBlock(  int l1, int m1, int l2, int m2,
                            int lmax, int mmax, int N, int order,
                            const std::vector<int>& Ms,
                            const std::vector<int>& mRows,
                            const BandedMatrix& expR,
                            std::function<maths::complex(int, int, int, int)> YlmYlm,
                            maths::Matrix m);

public:
    ExponentialPotential();

    double operator() (double x, double y, double z) const;
    void FillMatrix(maths::Matrix m, int N, const std::vector<int>& Ms = {0}, const std::vector<int>& mRows = {0});
    void FillMatrixGradX(maths::Matrix m, int N, int lmax, const std::vector<int>& Ms = {0}, const std::vector<int>& mRows = {0});
    void FillMatrixGradY(maths::Matrix m, int N, int lmax, const std::vector<int>& Ms = {0}, const std::vector<int>& mRows = {0});
    void FillMatrixGradZ(maths::Matrix m, int N, int lmax, const std::vector<int>& Ms = {0}, const std::vector<int>& mRows = {0});

    static std::string GetName();
    static tdse::Potential::Ptr_t Create(const nlohmann::json& observable);
    static bool Validate(const nlohmann::json& observable);
};