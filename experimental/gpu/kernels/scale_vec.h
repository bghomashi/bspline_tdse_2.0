#pragma once

#include "cl_gpu/cl_kernel.hpp"

auto scale_vec_src = StaticCode(
    __kernel void scale_vec(
        complex scale,                  // scale value
        __global complex* A,            // vector
        int i = get_global_id(0);     // each thread gets its own row
        
        A[i] = cmult(A[i], scale);
    }
);