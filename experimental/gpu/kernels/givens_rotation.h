#pragma once

#include "cl_gpu/cl_kernel.hpp"

// this is single threaded... I could not think of a better way :(
auto givens_rotation_col_src = StaticCode(
    __kernel void givens_rotation_col(
        __global complex* h,
        __global complex* cs,
        __global complex* sn,
        int h_cols,
        int k
    ) {
        int row1;
        int row2;
        // apply for ith column
        complex temp;
        // for each row in h
        for (int i = 0; i <= k-1; i++) {
            row1 = i*h_cols + k;
            row2 = (i+1)*h_cols + k;
            
            temp    =  cmult(cs[i],h[row1]) + cmult(sn[i],h[row2]);
            h[row2] = -cmult(sn[i],h[row1]) + cmult(cs[i],h[row2]);
            h[row1] = temp;
        }

        row1 = k*h_cols + k;
        row2 = (k+1)*h_cols + k;

        // update the next sin cos values for rotation
        double h1h1 = real(h[row1])*real(h[row1]) + imag(h[row1])*imag(h[row1]);    // |h1|^2
        double h2h2 = real(h[row2])*real(h[row2]);                                  // h2 should be real since it is the norm of a vector
        double t = sqrt(h1h1 + h2h2);                          
        complex cs_k = h[row1] / t;
        complex sn_k = h[row2] / t;                                                 // this is also real

        // eliminate H(i + 1, i)
        h[row1] = cmult(cs_k, h[row1]) + sn_k*h[row2];
        h[row2] = 0.0;

        cs[k] = cs_k;
        sn[k] = sn_k;
    }
);
