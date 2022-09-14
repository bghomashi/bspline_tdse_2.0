#pragma once

#include "cl_gpu/cl_kernel.hpp"

auto gram_schmidt_dot_partial_src = Kernel1D(
    complex,
    gram_schmidt_dot_partial, (
    __global const complex* Q,
    int a_col, int b_col, 
    int Qcols),
    int row = i;
    result = cmult((Q[row*Qcols + a_col]), conj(Q[row*Qcols + b_col]));
);

// Gram Schmidt
// assumes A is some array of complex scalars that a_index indexes into
// X is a matrix and we are orthogonalizing column (col_out)
//  against column (out)

auto gram_schmidt_iter_src = StaticCode(
    __kernel void gram_schmidt_iter(
        __global complex* h,                // scalar array
        __global complex* Q,                // vector to orthogonalize
        int a_index,                        // index into scalar
        int Qcols, int col_out, int col_in) {
            int i = get_global_id(0);       // element into vector
            Q[i*Qcols + col_out] -= cmult(h[a_index],Q[i*Qcols + col_in]);
    }
);