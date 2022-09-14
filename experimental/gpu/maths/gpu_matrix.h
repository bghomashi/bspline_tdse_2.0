// #ifndef __GPU_MATRIX_H__
// #define __GPU_MATRIX_H__

// #include "common/maths/matrix.h"
// #include "cl_gpu/cl_memory.hpp"
// // class PetscEngine;
// // class PetscEPS;
// // class PetscGMRES;

// class GPUMatrix : public maths::IMatrix {
//     // friend PetscEngine;
//     // friend PetscEPS;
//     // friend PetscGMRES;
// public:
//     typedef std::shared_ptr<GPUMatrix> Ptr_t;

//     Mat _petsc_mat;
//     int _row_start, _row_end;

//     // construct an empty matrix of size
//     GPUMatrix();
//     // construct a "banded" matrix of size rows*cols. Everything outside of 
//     // the last band is assumed to be 0
//     GPUMatrix(int rows, int cols, int numbands);
//     // copy contructor
//     // PetscMatrix(const PetscMatrix& o);
    
//     // destructor also deletes the matrix
//     ~GPUMatrix();

//     //  ------- basic operations described in "common/maths/matrix.h"  ------- 
//     maths::complex Get(int row, int col) const;
//     void Set(int row, int col, maths::complex val);
//     void Add(int row, int col, maths::complex val);
//     void Mult(const maths::Vector in, maths::Vector out);
//     void Scale(maths::complex factor);
//     void Duplicate(const maths::Matrix o);               // allocate AND copy values
//     void Zero();
//     void Copy(const maths::Matrix o);      

//     void FillBandedBlock(int bandwidth, int blocksize, int blockRow, int blockCol, maths::FuncOfRowCol element);
//     void FillBandedBlock(int bandwidth, int blocksize, maths::FuncOfRowCol element);
//     void FillBandedBlock(int bandwidth, int blockSize, int blockColOffset, maths::FuncOfRowCol element);

//     bool IsSymmetric(double tol) const;
//     bool IsAntiSymmetric(double tol) const;

//     void HermitianTranspose(maths::Matrix& transpose);
//     void Transpose(maths::Matrix& transpose);

//     void AYPX(maths::complex a, const maths::Matrix X);
//     void AXPY(maths::complex a, const maths::Matrix X);

//     // assemble the matrix
//     void AssembleBegin();
//     void AssembleEnd();
    
// };

// #endif