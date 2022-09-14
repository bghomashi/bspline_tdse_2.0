#pragma once

#include "cl_gpu/cl_kernel.hpp"

auto arnoldi_step_src = StaticCode(
    __kernel void arnoldi_step(
        __global const complex* A,      // matrix
        __global complex* Q,            // matrix of vectors
        int Qcols,                      // columns in Q
        int k) {                        // which column vector in Q to multiply
        int row = get_global_id(0);     // each thread gets its own row
        int Acols = get_global_size(0); // A is square
        // send output into next column of Q,
        // Q(k+1) = A*Q(k)
        for (int col = 0; col < Acols; col++)
            Q[row*Qcols + k+1] += cmult(A[row*Acols + col], Q[col*Qcols + k]);
    }
);