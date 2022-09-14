#pragma once

#include "cl_gpu/cl_kernel.hpp"

auto matvec_src = Kernel1D(
    complex,
    matmulvec,
    (__global const complex* A,     // matrix
    __global const complex* x,      // vector
    int Acols, int Xrows),
    int row = i;
    result = 0;
    for (int j = 0; j < Xrows; j++)
        result += cmult(A[row*Acols + j], x[j]);
);