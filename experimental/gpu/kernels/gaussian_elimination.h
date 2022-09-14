#pragma once

#include "cl_gpu/cl_kernel.hpp"

auto row_echelon_complex_src = Kernel1D(
    complex,
    row_echelon,
    (__global const complex* A, __global const complex* b),

    const int row = i;      // i - index of this kernel
    const int N = size;     // size - number of rows
    // copy b_i/A_ii into x_i for gaussian elimination step
    result = cdiv(b[row], A[row*N + row]);
);
auto gaussian_elimination_complex_src = Kernel1D(
    complex,
    gaussian_elimination,
    (__global const complex* A, const int col),

    const int row = i;      // i - index of this kernel
    const int N = size;     // size - number of rows
    const complex* x = ret_array;  // ret_array - values to be returned
    
    if (i < col)
        result -= cmult(x[col], cdiv(A[row*N + col],A[row*N + row]));
    else
        return;     // do nothing
);