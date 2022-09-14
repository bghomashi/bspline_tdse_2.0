#pragma once

#include "cl_gpu/cl_kernel.hpp"

auto dot_partial_src = Kernel1D(
    complex,
    dot_partial,
    (__global const complex* a, 
    __global const complex* b),
    result = a[i]*conj(b[i]);
);
