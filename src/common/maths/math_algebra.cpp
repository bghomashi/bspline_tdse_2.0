#include "common/maths/math_algebra.h"
#include "common/maths/math_common.h"
#include "common/maths/matrix.h"
#include "common/maths/vector.h"

namespace maths {
    void Mult(const Matrix M, const Vector in, Vector out) {
        M->Mult(in, out);
    }
    void Dot(const Vector a, const Vector b, complex& value) {
        a->Dot(b, value);
    }
    void AYPX(Matrix Y, complex a, const Matrix X) {
        Y->AYPX(a, X);
    }
    void AXPY(Matrix Y, complex a, const Matrix X) {
        Y->AXPY(a, X);
    }
    void AYPX(Vector Y, complex a, const Vector X) {
        Y->AYPX(a, X);
    }
    void AXPY(Vector Y, complex a, const Vector X) {
        Y->AXPY(a, X);
    }
}