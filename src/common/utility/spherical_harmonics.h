#pragma once

#include "common/maths/math_common.h"

maths::complex YlmXYlm(int l1, int m1, int l2, int m2);
maths::complex YlmYYlm(int l1, int m1, int l2, int m2);
maths::complex YlmZYlm(int l1, int m1, int l2, int m2);

maths::complex Ylm(int l, int m, double theta, double phi);
// *includes* complex conjugate of Y_l1^m1
maths::complex TripleIntegral(int l1, int l, int l2, int m1, int m, int m2);
