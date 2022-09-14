#pragma once

#include "cl_gpu/cl_kernel.hpp"

auto update_residual_src = StaticCode(
    __kernel void update_residual(
        __global complex* cs, 
        __global complex* sn, 
        __global complex* beta, 
        int k) {
        beta[k + 1] = -cmult(sn[k], beta[k]);
        beta[k]     =  cmult(cs[k], beta[k]);
        // error       = abs(beta[k+1]) / b_norm;
    }
);