#include "cl_kernels.h"

#include "complex_maths.h"
#include "zero_vec.h"
#include "sum_vector.h"
#include "dot_partial.h"
#include "gram_schmidt.h"
#include "matmul.h"
#include "arnoldi_step.h"
#include "scale_mat_col.h"
#include "gaussian_elimination.h"
#include "givens_rotation.h"
#include "update_residual.h"

void MathKernels::Initialize() {
    vector_program = CL::Program(
        complex_maths + 
        zero_vec_src + 
        // add_vec_src + 
        // scale_vec_src + 
        // AXPY_src + 
        // AYPX_src + 
        sum_src + 
        dot_partial_src
    );
    matrix_col_program = CL::Program(
        complex_maths + 
        gram_schmidt_dot_partial_src + 
        gram_schmidt_iter_src +
        arnoldi_step_src + 
        givens_rotation_col_src +
        scale_mat_col_src);

    matrix_vec_program = CL::Program(
        complex_maths + 
        update_residual_src + 
        matvec_src);
    
    gaussian_elimination_program = CL::Program(
        complex_maths + row_echelon_complex_src + gaussian_elimination_complex_src);


    sum = vector_program.CreateKernel("sum", 0);
    gram_schmidt_dot_partial    = matrix_col_program.CreateKernel("gram_schmidt_dot_partial", 1);
    gram_schmidt_iter           = matrix_col_program.CreateKernel("gram_schmidt_iter", 0);
    arnoldi_step                = matrix_col_program.CreateKernel("arnoldi_step", 0);
    scale_mat_col               = matrix_col_program.CreateKernel("scale_mat_col", 0);
    matmulvec                   = matrix_vec_program.CreateKernel("matmulvec", 1);
    givens_rotation_col         = matrix_col_program.CreateKernel("givens_rotation_col", 0);
    update_residual             = matrix_vec_program.CreateKernel("update_residual", 0);
    re_kernel                   = gaussian_elimination_program.CreateKernel("row_echelon", 1);
    ge_kernel                   = gaussian_elimination_program.CreateKernel("gaussian_elimination", 1);
}

void MathKernels::Destroy() {
    sum.~Kernel();
    gram_schmidt_iter.~Kernel();
    gram_schmidt_dot_partial.~Kernel();


    vector_program.~Program();
    matrix_col_program.~Program();
}










CL::Program MathKernels::vector_program;
CL::Program MathKernels::matrix_col_program;
CL::Program MathKernels::matrix_vec_program;
CL::Program MathKernels::gaussian_elimination_program;

CL::Kernel MathKernels::sum;
CL::Kernel MathKernels::gram_schmidt_iter, MathKernels::gram_schmidt_dot_partial, MathKernels::arnoldi_step, MathKernels::scale_mat_col;
CL::Kernel MathKernels::matmulvec;
CL::Kernel MathKernels::re_kernel, MathKernels::ge_kernel;
CL::Kernel MathKernels::givens_rotation_col, MathKernels::update_residual;