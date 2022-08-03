#ifndef __MATH_ALGEBRA_H__
#define __MATH_ALGEBRA_H__

#include "common/maths/math_common.h"

namespace maths {
    void Mult(const Matrix M, const Vector in, Vector out);
    void Dot(const Vector a, const Vector b, complex& value);
    void AYPX(Matrix Y, complex a, const Matrix X);
    void AXPY(Matrix Y, complex a, const Matrix X);
    void AYPX(Vector Y, complex a, const Vector X);
    void AXPY(Vector Y, complex a, const Vector X);
}

#endif