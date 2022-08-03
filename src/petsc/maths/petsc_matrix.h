#ifndef __PETSC_MATRIX_H__
#define __PETSC_MATRIX_H__

#include "common/maths/matrix.h"
#include <petsc.h>

class PetscEngine;
class PetscEPS;
class PetscGMRES;

class PetscMatrix : public maths::IMatrix {
    friend PetscEngine;
    friend PetscEPS;
    friend PetscGMRES;
public:
    typedef std::shared_ptr<PetscMatrix> Ptr_t;

    Mat _petsc_mat;
    int _row_start, _row_end;

    // construct an empty matrix of size
    PetscMatrix();
    // construct a "banded" matrix of size rows*cols. Everything outside of 
    // the last band is assumed to be 0
    PetscMatrix(int rows, int cols, int numbands);
    // copy contructor
    // PetscMatrix(const PetscMatrix& o);
    
    // destructor also deletes the matrix
    ~PetscMatrix();

    //  ------- basic operations described in "common/maths/matrix.h"  ------- 
    maths::complex Get(int row, int col) const;
    void Set(int row, int col, maths::complex val);
    void Add(int row, int col, maths::complex val);
    void Mult(const maths::Vector in, maths::Vector out);
    void Scale(maths::complex factor);
    void Duplicate(const maths::Matrix o);               // allocate AND copy values
    void Zero();
    void Copy(const maths::Matrix o);      

    void FillBandedBlock(int bandwidth, int blocksize, int blockRow, int blockCol, maths::FuncOfRowCol element);
    void FillBandedBlock(int bandwidth, int blocksize, maths::FuncOfRowCol element);
    void FillBandedBlock(int bandwidth, int blockSize, int blockColOffset, maths::FuncOfRowCol element);

    bool IsSymmetric(double tol) const;
    bool IsAntiSymmetric(double tol) const;

    void HermitianTranspose(maths::Matrix& transpose);
    void Transpose(maths::Matrix& transpose);

    void AYPX(maths::complex a, const maths::Matrix X);
    void AXPY(maths::complex a, const maths::Matrix X);

    // assemble the matrix
    void AssembleBegin();
    void AssembleEnd();

    // ------- petsc specific --------
    void Set(const std::vector<int>& rows, const std::vector<int>& cols, const maths::complex* value);
    int RowStart() const;
    int RowEnd() const;
        
};

#endif