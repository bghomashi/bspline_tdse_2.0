#pragma once

#include "cl_gpu/cl_kernel.hpp"


// C-style complex numbers
// - no operator overloading :(
static std::string complex_maths = std::string(StaticCode(
typedef double2 complex;

inline float  real(const complex a){
     return a.x;
}
inline float  imag(const complex a){
     return a.y;
}

inline float cabs(const complex a) {
    return (sqrt(a.x*a.x + a.y*a.y));
}

inline float carg(const complex a) {
    if(a.x > 0){
        return atan(a.y / a.x);

    }else if(a.x < 0 && a.y >= 0){
        return atan(a.y / a.x) + M_PI;

    }else if(a.x < 0 && a.y < 0){
        return atan(a.y / a.x) - M_PI;

    }else if(a.x == 0 && a.y > 0){
        return M_PI/2.;

    }else if(a.x == 0 && a.y < 0){
        return -M_PI/2.;

    }else{
        return 0;
    }
}

inline complex conj(const complex a) {
    return (complex)(a.x, -a.y);
}

inline complex cmult(const complex a, const complex b) {
    return (complex)(a.x*b.x - a.y*b.y, a.x*b.y + a.y*b.x);
}

inline complex cdiv(const complex a, const complex b) {
    return (complex)((a.x*b.x + a.y*b.y)/(b.x*b.x + b.y*b.y), (a.y*b.x - a.x*b.y)/(b.x*b.x + b.y*b.y));
}

inline complex csqrt(const complex a){
    return (complex)( sqrt(cabs(a)) * cos(carg(a)/2.),  sqrt(cabs(a)) * sin(carg(a)/2.));
}
)) + "#define I ((complex)(0.0, 1.0))\n\n";