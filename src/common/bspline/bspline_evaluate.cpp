#include "common/bspline/bspline.h"

#include <complex>

using namespace bspline;
std::vector<maths::complex> BSpline::getBSpline(const std::vector<double>& x, int i, int dn) const {
    std::vector<maths::complex> r(x.size());

    for (int j = 0; j < (int)x.size(); j++) {
        r[j] = bspline(x[j], i, dn);
    }
    return r;
}

maths::complex BSpline::bspline(double x, int bs, int dn) const {
    if (bs < 0 || bs >= _numBSplines) return 0;

    int i = whichInterval(x);						// which grid interval is 'x'
    int interval = i - bs + _order-1;					// which interval is this in the bspline

    if (interval < 0 || interval > _order-1) return 0;


    maths::complex sum = 0;
    double xx = 1;
    double invInt = 1. / (_grid[i + 1] - _grid[i]);
    double r = (x - _grid[i]) * invInt;

    const auto& coeff = _bsCoeffs[bs];

    for (int d = dn; d < _order; d++) {
        // std::cout << "interval: " << interval << " degree: " << d << " = " << coeff(interval, d) << std::endl;
        sum += _partialFactorial[d + dn * DIMFACT] * coeff(interval, d) * xx;
        xx *= r;
    }
    // 
    double invIntN = 1.;
    for (int k = 0; k < dn; k++)
        invIntN *= invInt;

    return sum * invIntN;
}
maths::complex BSpline::bspline(maths::complex x, int bs, int dn) const {
    if (bs < 0 || bs >= _numBSplines) return 0;

    int i = whichInterval(_ecs.x(x));	// which grid interval is 'x'
    int interval = i - bs + _order-1;					// which interval is this in the bspline

    if (interval < 0 || interval > _order-1) return 0;


    maths::complex sum = 0;
    maths::complex xx = 1;
    // can be optimized??
    maths::complex invInt = 1. / (_ecs_grid[i + 1] - _ecs_grid[i]);
    maths::complex r = (x - _ecs_grid[i]) * invInt;		// << this should be real


    const auto& coeff = _bsCoeffs[bs];

    for (int d = dn; d < _order; d++) {
        // std::cout << "interval: " << interval << " degree: " << d << " = " << coeff(interval, d) << std::endl;
        sum += _partialFactorial[d + dn*DIMFACT] * coeff(interval, d) * xx;
        xx *= r;
    }

    // 
    maths::complex invIntN = 1.;
    for (int k = 0; k < dn; k++)
        invIntN *= invInt;

    return sum * invIntN;
}


// - functions for a single vector
// std::vector<complex> Basis::FunctionEvaluate(const std::vector<double>& x, const std::vector<double>& fc, int dn) const {
// 	std::vector<complex> r(x.size());

// 	for (int j = 0; j < (int)x.size(); j++) {
// 		r[j] = FunctionEvaluate(x[j], fc, dn);
// 	}
// 	return r;
// }

maths::complex BSpline::FunctionEvaluate(double x, const std::vector<maths::complex>& fc, int dn) const {
    int Bsmi, Bsma;
    int nskip;

    maths::complex total = 0.0;
    int i = whichInterval(x);						// where is 'x'

    if (i < 0 || i > _nodes - 2) return 0.;			// if it is outside the grid return 0

    Bsmi = 0;
    nskip = 0;

    if (_skipFirst) {
        Bsmi = 1;
        nskip = 1;
    }

    Bsmi = std::max(Bsmi, i);
    Bsma = std::min(_numBSplines, i + _order);

    for (int Bs = Bsmi; Bs < Bsma; Bs++) {
        if (Bs - nskip >= fc.size())                   // BUGFIX:: I think this ONLY happens if we are skipping the LAST spline??
            continue;
        total += bspline(x, Bs, dn) * fc[Bs - nskip];
    }
    return total;
}



// - functions for a single vector
std::vector<maths::complex> BSpline::FunctionEvaluate(const std::vector<double>& x, const std::vector<maths::complex>& fc, int dn) const {
    std::vector<maths::complex> r(x.size());

    for (int j = 0; j < (int)x.size(); j++) {
        r[j] = FunctionEvaluate(x[j], fc, dn);
    }
    return r;
}
