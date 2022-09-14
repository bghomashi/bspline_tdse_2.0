#ifndef __GPU_SOLVER_H__
#define __GPU_SOLVER_H__

#include "common/maths/math_common.h"
#include "common/maths/gmres_solver.h"
#include "cl_gpu/cl_memory.hpp"

class GPUSolver : public maths::IGMRESSolver {
    typedef CL::Vector1D<maths::complex> Vector1D;
    
    class SqrMatrix {
    public:
        CL::Vector1D<maths::complex> elements;
        int rows, cols;

        SqrMatrix() {
            rows = 0; 
            cols = 0;
        }
        SqrMatrix(int rows, int cols) : rows(rows), cols(cols) {
            elements.resize(rows*cols);
        }
        void Resize(int rows, int cols) {
            elements.resize(rows*cols);
        }

        maths::complex& operator() (int r, int c) {
            return elements(r*cols + c);
        }
    };



    int restart_iter, max_iter;

    Vector1D gs_intermediate_vector;
    Vector1D beta, cs, sn;
    SqrMatrix h, Q;
public:

    GPUSolver(int restart_iter = 500, int max_iter = 10000);
    ~GPUSolver();

    void SetBlockedPC(int blocks);
    bool Solve(const maths::Matrix A, const maths::Vector b, maths::Vector x);



    void ArnoldiStep(Matrix& A, Matrix& Q, int k);
    void GramSchmidt(Matrix& h, Matrix& Q, int k);
    void ApplyGivensRotation(Matrix& h, Vector& cs, Vector& sn, int k);
    void UpdateResidual(Vector& beta, Vector& cs, Vector& sn, int k);
};

#endif
