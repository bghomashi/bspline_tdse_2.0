#include "common/bspline/bspline.h"
#include "common/bspline/bspline_common.h"

#include <algorithm>
#include <cmath>
#include <iostream>

using namespace bspline;
bool GaussQuadrature::Initialize() {
    // initialize gaussian points/weigths
    _points.resize(NGAUSS);
    _weights.resize(NGAUSS);

    double p1, p2, p3, pp, z, z1;

    for (int i = 1; i <= std::floor((NGAUSS + 1) / 2.); i++) {
        z = std::cos(maths::Pi * (i - .25) / (NGAUSS + .5));
        do {
            p1 = 1.0;
            p2 = 0.0;
            for (int j = 1; j <= NGAUSS; j++) {
                p3 = p2;
                p2 = p1;
                p1 = ((2.0 * j - 1.0) * z * p2 - (j - 1.0) * p3) / double(j);
            }
            pp = double(NGAUSS) * (z * p1 - p2) / (z * z - 1.0);
            z1 = z;
            z = z1 - p1 / pp;
        } while (std::abs(z - z1) > GAUSS_POINTS_CONVERGENCE_THRESHOLD);

        _points[i - 1] = (1. - z) / 2.;
        _points[NGAUSS - i] = (1. + z) / 2.;
        _weights[i - 1] = 1. / ((1. - z * z) * pp * pp);
        _weights[NGAUSS - i] = 1. / ((1. - z * z) * pp * pp);
    }

    return true;
}


maths::complex GaussQuadrature::Integrate(maths::complex xmin, maths::complex xmax, std::function<maths::complex(maths::complex)> f) const {
    maths::complex width = xmax - xmin;
    maths::complex sum = 0.;
    maths::complex c = 0.;
    maths::complex t = 0., y = 0., elem = 0.;
    
    if (std::abs(width) < NOD_THRESHOLD) return 0;	// spacing is too small, return 0;
    for (int i = 0; i < NGAUSS; i++) {
        elem = f(xmin + width * _points[i]) * _weights[i];
        y = elem - c;
        t = sum + y;
        c = (t-sum) - y;
        sum = t;
    }


    return sum*width;
}