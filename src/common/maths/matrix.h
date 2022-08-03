#ifndef __MATRIX_H__
#define __MATRIX_H__

#include "common/maths/math_common.h"

// Establishes the interface of any Matrix implementation to 
// be used by the TDSE/TISE methods and/or to compute any
// observables.

namespace maths {
    class IMatrix {
    protected:
        int _rows, _cols;
    public:
        // retrieve the rows and cols of this matrix
        int Rows() const {
            return _rows;
        }
        int Cols() const {
            return _cols;
        }

        // get the value at some row, col
        virtual complex Get(int row, int col) const = 0;
        // set the value at some row, col
        virtual void Set(int row, int col, complex val) = 0;
        // add a value to the element at some row, col
        virtual void Add(int row, int col, complex val) = 0;
        // multiply this matrix on to the vector in. The result
        // is stored in the vector out
        virtual void Mult(const Vector in, Vector out) = 0;
        // scale all elements of the matrix
        virtual void Scale(complex factor) = 0;
        // make matrix have the same struct AND copy values
        // from matrix o.
        virtual void Duplicate(const Matrix o) = 0;               
        // zero the entire matrix
        virtual void Zero() = 0;
        // copy the values from o. Does no allocate structure.
        virtual void Copy(const Matrix o) = 0;

        // test if a this matrix is symmetric to some tolerance
        virtual bool IsSymmetric(double tol = 1e-17) const = 0;
        // test if a this matrix is anti-symmetric to some tolerance
        virtual bool IsAntiSymmetric(double tol = 1e-17) const = 0;
        // get the hermitian transpose
        virtual void HermitianTranspose(maths::Matrix& transpose) = 0;
        // get the transpose
        virtual void Transpose(maths::Matrix& transpose) = 0;

        // multiply (matrix) this vector by (scalar) a then add (matrix) X
        // results replaces this
        virtual void AYPX(complex a, Matrix X) = 0;
        // multiply (matrix) X by (scalar) a then add (matrix) this
        // results replaces this
        virtual void AXPY(complex a, Matrix X) = 0;

        // The next three functions fill a banded block in a block matrix such as 
        // finite difference representation of poison equation. Each block has a 
        // bandwidth and a total size. 
        // FuncOfRowCol is a function taking the row and col as parameters and 
        // returns a complex value.

        // specify a particular block (blockRow, blockCol) to fill. The blockRow 
        // and blockCol are in units of the blockSize.
        virtual void FillBandedBlock(int bandwidth, int blocksize, int blockRow, int blockCol, FuncOfRowCol element) = 0;
        // fills main diagonal blocks in the matrix
        virtual void FillBandedBlock(int bandwidth, int blocksize, FuncOfRowCol element) = 0;
        // fills sub diagonal blocks in the matrix. Shifted by blockColOffset.
        virtual void FillBandedBlock(int bandwidth, int blockSize, int blockColOffset, FuncOfRowCol element) = 0;


        // for PETSc...
        virtual void AssembleBegin() {};
        virtual void AssembleEnd() {};
    };
}

#endif