#pragma once

#include "cl_gpu/cl_kernel.hpp"

static std::string ylmylm_maths = StaticCode(
    // evaluates:  <l1 m1| Xi/r |l2 m2> where Xi = x, y, or z
    double YlmXYlm(int l1, int m1, int l2, int m2) {
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
    double YlmYYlm(int l1, int m1, int l2, int m2) {                   // SIGN FLIP??
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
    double YlmZYlm(int l1, int m1, int l2, int m2) {
        if (m1 != m2) return 0;
        if (l1 == l2 + 1) 
            return sqrt( (l2-m2+1.)*(l2+m2+1.)/(2.*l2+1.)/(2.*l2+3.) );
        
        if (l1 == l2 - 1) 
            return sqrt( (l2-m2)*(l2+m2)/(2.*l2+1)/(2.*l2-1) );

        return 0.;
    }
);

