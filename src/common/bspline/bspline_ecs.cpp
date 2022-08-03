#include "common/bspline/bspline_ecs.h"
#include "common/bspline/bspline_common.h"
#include <cmath>

using namespace bspline;

maths::complex ECS::q(double x) const {
    if (x <= r0)
        return 1;
    return std::exp(I*theta);
}
maths::complex ECS::R(double x) const {
    if (x <= r0)
        return x;
    return r0 + (x-r0)*std::exp(I*theta);
}
double ECS::x(maths::complex R) const {
    if (std::imag(R) == 0)
        return std::real(R);
    return r0 + std::real((R - r0)*std::exp(-I*theta));
}
