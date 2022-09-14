#include "gpu/maths/gpu_gmres.h"
#include <iostream>
#include "gpu/kernels/cl_kernels.h"
#include "gpu/maths/gpu_vector.h"
#include "gpu/maths/gpu_common.h"

using namespace maths;


GPUSolver::GPUSolver(int restart_iter, int max_iter) : restart_iter(restart_iter), max_iter(max_iter) {
    beta.resize(restart_iter+1);
    cs.resize(restart_iter);
    sn.resize(restart_iter);
    h.Resize(restart_iter+1, restart_iter);
    Q.Resize(restart_iter, restart_iter+1);
}
GPUSolver::~GPUSolver() {
}

void GPUSolver::SetBlockedPC(int blocks) {
}

bool GPUSolver::Solve(const Matrix A, const Vector vb, Vector vx) {
    auto b = GPUCast(vb);
    complex norm = 0;
    for (int i = 0; i < N; i++)
        norm += b->_gpu_vec(i)*conj(b(i));


    norm = sqrt(norm);
    for (int i = 0; i < N; i++) {
        Q(i, 0) = b(i) / norm;
    }
    beta(0) = norm;

    
        ArnoldiStep(A, Q, k);
        GramSchmidt(h, Q, k);
        ApplyGivensRotation(h, cs, sn, k);
        UpdateResidual(beta, cs, sn, k);

    return true;
}






void GPUSolver::ArnoldiStep(Matrix& A, Matrix& Q, int k) {
    MathKernels::arnoldi_step.setArg(0, A.elements);
    MathKernels::arnoldi_step.setArg(1, Q.elements);
    MathKernels::arnoldi_step.setArg(2, Q.cols);
    MathKernels::arnoldi_step.setArg(3, k);
    MathKernels::arnoldi_step.Execute(A.rows);
}

void GPUSolver::GramSchmidt(Matrix& h, Matrix& Q, int k) {
    gs_intermediate_vector.resize(Q.rows);

    for (int i = 0; i <= k; i++) {         // past k vectors
        // do partial dot product
        gs_intermediate_vector = MathKernels::gram_schmidt_dot_partial(
            Q.elements,         // matrix vectors
            k+1,                // column to project
            i,                  // column to project on
            Q.cols);            // how many total columns in 'Q'
        
        // reduce
        int len = gs_intermediate_vector.size();
        while (len > 1) {
            int sums = len / 2;
            MathKernels::sum.setArg(0, gs_intermediate_vector._buffer);
            MathKernels::sum.setArg(1, len);
            MathKernels::sum.Execute(sums);
            len -= sums;
        }
        h.elements._buffer.Copy(
            gs_intermediate_vector._buffer, 0,      // where value of dot-product is stored
            (i*h.cols + k)*sizeof(complex),         // index to which scalar in 'h' in bytes
            sizeof(complex));                       // how big is element in bytes

        // // orthogonalize against the k'th colum
        MathKernels::gram_schmidt_iter.setArg(0, h.elements._buffer);   // array storing projections
        MathKernels::gram_schmidt_iter.setArg(1, Q.elements._buffer);   // matrix of vectors
        MathKernels::gram_schmidt_iter.setArg(2, i*h.cols + k);       // index to which scalar in 'h'
        MathKernels::gram_schmidt_iter.setArg(3, Q.cols);               // how many total columns in 'Q'
        MathKernels::gram_schmidt_iter.setArg(4, k+1);                  // column to orthogonalize
        MathKernels::gram_schmidt_iter.setArg(5, i);                    // column to orthogonalize against

        MathKernels::gram_schmidt_iter.Execute(Q.rows);  
    }

    // normalize Q(k+1)
    // do partial dot product
    gs_intermediate_vector = MathKernels::gram_schmidt_dot_partial(
        Q.elements,         // matrix vectors
        k+1,                // column to project
        k+1,                // column to project on
        Q.cols);            // how many total columns in 'Q'

    // reduce
    int len = gs_intermediate_vector.size();
    while (len > 1) {
        int sums = len / 2;
        MathKernels::sum.setArg(0, gs_intermediate_vector._buffer);
        MathKernels::sum.setArg(1, len);
        MathKernels::sum.Execute(sums);
        len -= sums;
    }

    // unfortunate to read back - probably not worth writing a kernel just for this?
    // read back - sqrt - write to h
    complex L;
    gs_intermediate_vector._buffer.Read(&L, 0, sizeof(complex));
    L = sqrt(L);
    assert(L.imag() < 1e-10 && "imaginary part of 2-norm is nonzero.");
    L.imag(0);                  // norm of a vect is real
    h.elements._buffer.Write(&L, ((k+1)*h.cols + k)*sizeof(complex), sizeof(complex));
    
    if (std::abs(L) > ESPILON) {
        MathKernels::scale_mat_col.setArg(0, 1./L);
        MathKernels::scale_mat_col.setArg(1, Q.elements);
        MathKernels::scale_mat_col.setArg(2, Q.cols);
        MathKernels::scale_mat_col.setArg(3, k+1);
        MathKernels::scale_mat_col.Execute(Q.rows);
    }
}

void GPUSolver::ApplyGivensRotation(Matrix& h, Vector& cs, Vector& sn, int k) {
    MathKernels::givens_rotation_col.setArg(0, h.elements);
    MathKernels::givens_rotation_col.setArg(1, cs.elements);
    MathKernels::givens_rotation_col.setArg(2, sn.elements);
    MathKernels::givens_rotation_col.setArg(3, h.cols);
    MathKernels::givens_rotation_col.setArg(4, k);
    MathKernels::givens_rotation_col.Execute(1);
}

void GPUSolver::UpdateResidual(Vector& beta, Vector& cs, Vector& sn, int k) {
    MathKernels::update_residual.setArg(0, cs.elements);
    MathKernels::update_residual.setArg(1, sn.elements);
    MathKernels::update_residual.setArg(2, beta.elements);
    MathKernels::update_residual.setArg(3, k);
    MathKernels::update_residual.Execute(1);
}
