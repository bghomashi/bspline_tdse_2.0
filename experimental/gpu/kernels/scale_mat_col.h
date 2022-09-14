#pragma once

#include "cl_gpu/cl_kernel.hpp"

auto scale_mat_col_src = StaticCode(
    __kernel void scale_mat_col(
        complex scale,                  // scale value
        __global complex* A,            // matrix
        int Acols,                      // columns in A
        int k) {                        // which column vector in Q to multiply
        int row = get_global_id(0);     // each thread gets its own row
        
        A[row*Acols + k] = cmult(A[row*Acols + k], scale);
    }
);