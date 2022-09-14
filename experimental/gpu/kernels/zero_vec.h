#pragma once

#include "cl_gpu/cl_kernel.hpp"

auto zero_vec_src = Kernel1D(
    complex,
    zero_vec,
    (int a),
    result = (complex)(0., 0.);
);