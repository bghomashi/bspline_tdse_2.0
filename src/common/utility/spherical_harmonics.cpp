#include "common/utility/spherical_harmonics.h"
#include <cmath>
#include <complex>
#include "common/utility/wigner3j/wigner3j.h"

using namespace std::complex_literals;

maths::complex YlmXYlm(int l1, int m1, int l2, int m2) {
    if (m1 == m2 + 1 && l1 == l2 + 1)
        return -0.5*sqrt( (l2+m2+1.)*(l2+m2+2.)/(2.*l2+1.)/(2.*l2+3.) );

    if (m1 == m2 + 1 && l1 == l2 - 1)
        return 0.5*sqrt( (l2-m2)*(l2-m2-1.)/(2.*l2-1.)/(2.*l2+1.) );

    if (m1 == m2 - 1 && l1 == l2 + 1)
        return 0.5*sqrt( (l2-m2+1.)*(l2-m2+2.)/(2.*l2+1.)/(2.*l2+3.) );

    if (m1 == m2 - 1 && l1 == l2 - 1)
        return -0.5*sqrt( (l2+m2)*(l2+m2-1.)/(2.*l2-1.)/(2.*l2+1.) );

    return 0.;
}
maths::complex YlmYYlm(int l1, int m1, int l2, int m2) {                   // SIGN FLIP??
    if (m1 == m2 + 1 && l1 == l2 + 1)
        return 0.5i*sqrt( (l2+m2+1.)*(l2+m2+2.)/(2.*l2+1.)/(2.*l2+3.) );

    if (m1 == m2 + 1 && l1 == l2 - 1)
        return -0.5i*sqrt( (l2-m2)*(l2-m2-1.)/(2.*l2-1.)/(2.*l2+1.) );

    if (m1 == m2 - 1 && l1 == l2 + 1)
        return 0.5i*sqrt( (l2-m2+1.)*(l2-m2+2.)/(2.*l2+1.)/(2.*l2+3.) );

    if (m1 == m2 - 1 && l1 == l2 - 1)
        return -0.5i*sqrt( (l2+m2)*(l2+m2-1.)/(2.*l2-1.)/(2.*l2+1.) );

    return 0.;
}
maths::complex YlmZYlm(int l1, int m1, int l2, int m2) {
    if (m1 != m2) return 0;
    if (l1 == l2 + 1) 
        return sqrt( (l2-m2+1.)*(l2+m2+1.)/(2.*l2+1.)/(2.*l2+3.) );
    
    if (l1 == l2 - 1) 
        return sqrt( (l2-m2)*(l2+m2)/(2.*l2+1)/(2.*l2-1) );
    

    return 0.;
}




maths::complex Ylm(int l, int m, double theta, double phi) {
    return std::sph_legendre(l, m, theta)*std::exp(1i*phi*(double)m);
}

maths::complex TripleIntegral(int l1, int l, int l2, int m1, int m, int m2) {
	double sign = (m1 % 2 == 0 ? 1. : -1.);
	return sign * 
		std::sqrt((2.*l1 + 1.)*(2.*l + 1.)*(2.*l2 + 1.)/(4.*maths::Pi)) *
		Wigner3j::ThreeJd(l1, l, l2, 0, 0, 0) * 
		Wigner3j::ThreeJd(l1, l, l2, -m1, m, m2);
}