#include "common/bspline/bspline.h"



static bspline::BSpline s_basis;

namespace bspline {
    bool Basis::Validate(const nlohmann::json& basis) {
        return s_basis.Validate(basis);
    }
    void Basis::Load(const nlohmann::json& basis, bool ecs_off) {
        s_basis.Load(basis, ecs_off);
    }
    maths::complex Basis::bspline(maths::complex x, int bs, int dn) {
        return s_basis.bspline(x, bs, dn);
    }
    maths::complex Basis::bspline(double x, int bs, int dn) {
        return s_basis.bspline(x, bs, dn);
    }
    std::vector<maths::complex> Basis::GetBSpline(const std::vector<double>& x, int bs, int dn) {
        return s_basis.getBSpline(x, bs, dn);
    }

    maths::complex Basis::FunctionEvaluate(double x, const std::vector<maths::complex>& fc, int dn) {
        return s_basis.FunctionEvaluate(x, fc, dn);
    }
    std::vector<maths::complex> Basis::FunctionEvaluate(const std::vector<double>& x, const std::vector<maths::complex>& fc, int dn) {
        return s_basis.FunctionEvaluate(x, fc, dn);
    }

    maths::complex Basis::Integrate(int bs1, int bs2, int dn1, int dn2) {
        return s_basis.Integrate(bs1, bs2, dn1, dn2);
    }
    maths::complex Basis::Integrate(int bs1, int bs2, std::function<maths::complex(maths::complex)> f, int dn1, int dn2) {
        return s_basis.Integrate(bs1, bs2, f, dn1, dn2);
    }
    maths::complex Basis::Integrate(double xmin, double xmax, int bs1, int bs2, std::function<maths::complex(maths::complex)> f, int dn1, int dn2) {
        return s_basis.Integrate(xmin, xmax, bs1, bs2, f, dn1, dn2);
    }

    const std::vector<double>& Basis::GetGrid() {
        return s_basis.getGrid();
    }
    int Basis::GetNumBSplines()  {
        return s_basis.getNumBSplines();
    }
    int Basis::GetOrder()  {
        return s_basis.getOrder();
    }
    int Basis::GetNumNodes() {
        return s_basis.getNumNodes();
    }
    double Basis::GetECS_R0() {
        return s_basis.getECS_R0();
    }
    double Basis::GetECS_Theta() {
        return s_basis.getECS_Theta();
    }
    double Basis::GetXmin() {
        return s_basis.getXmin();
    }
    double Basis::GetXmax() {
        return s_basis.getXmax();
    }
}