#pragma once

#include "cl_gpu/cl_kernel.hpp"

class MathKernels {
public:
    static CL::Program vector_program;
    static CL::Program matrix_col_program;
    static CL::Program matrix_vec_program;
    static CL::Program gaussian_elimination_program;
    
    static CL::Kernel sum;
    static CL::Kernel gram_schmidt_iter, gram_schmidt_dot_partial, arnoldi_step, scale_mat_col;
    static CL::Kernel matmulvec;
    static CL::Kernel re_kernel, ge_kernel;
    static CL::Kernel givens_rotation_col, update_residual;


    static void Initialize();
    static void Destroy();
};