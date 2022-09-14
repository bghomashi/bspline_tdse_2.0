// #include "petsc/maths/petsc_matrix.h"
// #include "petsc/maths/petsc_common.h"
// #include <petscmat.h>

// using namespace maths;

// PetscMatrix::PetscMatrix() {
//     _rows = 0; _cols = 0;
//     _row_start = 0; _row_end = 0;
//     _petsc_mat = 0;
// }


// PetscMatrix::PetscMatrix(int rows, int cols, int numbands) {
//     PetscErrorCode ierr;
//     ierr = MatCreate(PETSC_COMM_WORLD,&_petsc_mat);PETSCASSERT(ierr);
//     ierr = MatSetSizes(_petsc_mat,PETSC_DECIDE,PETSC_DECIDE,rows,cols);PETSCASSERT(ierr);
//     ierr = MatMPIAIJSetPreallocation(_petsc_mat, numbands, NULL, numbands, NULL);PETSCASSERT(ierr);
//     ierr = MatSetOption(_petsc_mat,MAT_IGNORE_ZERO_ENTRIES,PETSC_TRUE);PETSCASSERT(ierr);
//     ierr = MatSetOption(_petsc_mat,MAT_NEW_NONZERO_LOCATIONS,PETSC_FALSE);PETSCASSERT(ierr);
//     ierr = MatSetOption(_petsc_mat,	MAT_NEW_NONZERO_LOCATION_ERR,PETSC_TRUE);PETSCASSERT(ierr);
//     ierr = MatSetFromOptions(_petsc_mat);PETSCASSERT(ierr);
//     ierr = MatSetUp(_petsc_mat);PETSCASSERT(ierr);
//     ierr = MatGetOwnershipRange(_petsc_mat,&_row_start,&_row_end); PETSCASSERT(ierr);
//     _rows = rows; _cols = cols;
// }
        
// PetscMatrix::~PetscMatrix() {
//     MatDestroy(&_petsc_mat);
//     _petsc_mat = 0;
// }

// void PetscMatrix::AssembleBegin() {
//     PetscErrorCode ierr = MatAssemblyBegin(_petsc_mat,MAT_FINAL_ASSEMBLY);PETSCASSERT(ierr);
// }
// void PetscMatrix::AssembleEnd() {
//     PetscErrorCode ierr = MatAssemblyEnd(_petsc_mat,MAT_FINAL_ASSEMBLY);PETSCASSERT(ierr);
// }


// complex PetscMatrix::Get(int row, int col) const {
//     return 0.;              // no obvious way to do this in PETSc?
// }
// void PetscMatrix::Set(int row, int col, complex val) {
//     PetscErrorCode ierr = MatSetValue(_petsc_mat, row, col, val, INSERT_VALUES); PETSCASSERT(ierr);
// }
// void PetscMatrix::Add(int row, int col, complex val) {
//     PetscErrorCode ierr = MatSetValue(_petsc_mat, row, col, val, ADD_VALUES); PETSCASSERT(ierr);
// }
// void PetscMatrix::Mult(const Vector in, Vector out) {
//     auto a = PetscCast(in);
//     auto b = PetscCast(out);

//     PetscErrorCode ierr = MatMult(_petsc_mat, a->_petsc_vec, b->_petsc_vec); PETSCASSERT(ierr);
// }
// void PetscMatrix::Scale(complex factor) {
//     PetscErrorCode ierr = MatScale(_petsc_mat, factor); PETSCASSERT(ierr);
// }
// void PetscMatrix::Duplicate(const Matrix o) {
//     PetscErrorCode ierr;
//     // destroy previous matrix if any
//     MatDestroy(&_petsc_mat);
//     _petsc_mat = 0;

//     auto from = PetscCast(o);
//     ierr = MatDuplicate(from->_petsc_mat, MAT_COPY_VALUES, &_petsc_mat); PETSCASSERT(ierr);
//     ierr = MatGetOwnershipRange(_petsc_mat,&_row_start,&_row_end); PETSCASSERT(ierr);
//     _rows = from->_rows; _cols = from->_cols;
// }
// void PetscMatrix::Zero() {
//     MatZeroEntries(_petsc_mat);
// }

// void PetscMatrix::Copy(const Matrix o) {
//     auto from = std::dynamic_pointer_cast<PetscMatrix>(o);
//     MatCopy(from->_petsc_mat, _petsc_mat, SAME_NONZERO_PATTERN);
// }

// void PetscMatrix::Set(const std::vector<int>& rows, const std::vector<int>& cols, const complex* values) {
//     MatSetValues(_petsc_mat, rows.size(), rows.data(), cols.size(), cols.data(), values, INSERT_VALUES);
// }

// void PetscMatrix::FillBandedBlock(int bandwidth, int blocksize, int blockRow, int blockCol, std::function<complex(int, int)> element) {
//     std::vector<int> cols; cols.reserve(2*bandwidth+1);
//     std::vector<complex> values; values.reserve(2*bandwidth+1);

//     if  (blockRow*blocksize >= _row_end || 
//         (blockRow+1)*blocksize <= _row_start)
//         return;                         // outside of our range

//     int start = std::max(_row_start, blockRow*blocksize);
//     int end = std::min(_row_end, (blockRow+1)*blocksize);


//     for (int r = start; r < end; r++) {
//         if (blockCol > (_cols/blocksize) ||
//             blockCol < 0)
//             continue;

//         int colStart = std::max(blockCol*blocksize, blockCol*blocksize+(r-blockRow*blocksize)-bandwidth), 
//             colEnd = std::min((blockCol+1)*blocksize-1, blockCol*blocksize+(r-blockRow*blocksize)+bandwidth);

//         for (int c = colStart; c <= colEnd; c++) {
//             cols.push_back(c);
//             values.push_back(element(r, c));
//         }
//         Set({r}, cols, values.data());
//         cols.clear();
//         values.clear();
//     }
// }

// void PetscMatrix::FillBandedBlock(int bandwidth, int blocksize, std::function<complex(int, int)> element) {
//     std::vector<int> cols; cols.reserve(2*bandwidth+1);
//     std::vector<complex> values; values.reserve(2*bandwidth+1);

//     for (int r = _row_start; r < _row_end; r++) {
//         int block = std::floor(r / blocksize);

//         int colStart = std::max(block*blocksize, r-bandwidth), 
//             colEnd = std::min((block+1)*blocksize-1, r+bandwidth);

//         for (int c = colStart; c <= colEnd; c++) {
//             cols.push_back(c);
//             values.push_back(element(r, c));
//         }
//         Set({r}, cols, values.data());
//         cols.clear();
//         values.clear();
//     }

//     AssembleBegin();
//     AssembleEnd();
// }
// void PetscMatrix::FillBandedBlock(int bandwidth, int blockSize, int blockColOffset, std::function<complex(int, int)> element) {
//     std::vector<int> cols; cols.reserve(2*bandwidth+1);
//     std::vector<complex> values; values.reserve(2*bandwidth+1);

//     for (int r = _row_start; r < _row_end; r++) {
//         int block = std::floor(r / blockSize);
//         if (block+blockColOffset > (_cols/blockSize) ||
//             block+blockColOffset < 0)
//             continue;

//         int colStart = std::max((block+blockColOffset)*blockSize, (r+blockColOffset*blockSize)-bandwidth), 
//             colEnd = std::min((block+blockColOffset+1)*blockSize-1, (r+blockColOffset*blockSize)+bandwidth);

//         if (colStart >= _cols) break;

//         for (int c = colStart; c <= colEnd; c++) {
//             cols.push_back(c);
//             values.push_back(element(r, c));
//         }
//         Set({r}, cols, values.data());
//         cols.clear();
//         values.clear();
//     }
    
//     AssembleBegin();
//     AssembleEnd();
// }
// bool PetscMatrix::IsSymmetric(double tol) const {
//     PetscBool result;
//     MatIsSymmetric(_petsc_mat, tol, &result);
//     return (result == PETSC_TRUE);
// }
// bool PetscMatrix::IsAntiSymmetric(double tol) const {
//     PetscReal norm;
//     Mat B;

//     MatTranspose(_petsc_mat, MAT_INITIAL_MATRIX, &B);
//     MatAXPY(B, 1.0, _petsc_mat, SAME_NONZERO_PATTERN);
//     MatNorm(B, NORM_FROBENIUS, &norm);


//     MatDestroy(&B);

//     if (norm <= tol)
//         return true;

//     // std::cout << "2-norm=" << norm <<std::endl;
//     return false;
// }

// void PetscMatrix::HermitianTranspose(maths::Matrix& transpose) {
//     auto B = PetscCast(transpose);
//     PetscErrorCode ierr = MatHermitianTranspose(_petsc_mat, MAT_REUSE_MATRIX, &B->_petsc_mat); PETSCASSERT(ierr);
// }
// void PetscMatrix::Transpose(maths::Matrix& transpose) {
//     auto B = PetscCast(transpose);
//     PetscErrorCode ierr = MatTranspose(_petsc_mat, MAT_REUSE_MATRIX, &B->_petsc_mat); PETSCASSERT(ierr);
// }
// void PetscMatrix::AYPX(maths::complex a, const Matrix X) {
//     auto y = this;
//     auto x = PetscCast(X);

//     MatAYPX(y->_petsc_mat, a, x->_petsc_mat, SUBSET_NONZERO_PATTERN);
// }
// void PetscMatrix::AXPY(maths::complex a, const Matrix X) {
//     auto y = this;
//     auto x = PetscCast(X);

//     MatAXPY(y->_petsc_mat, a, x->_petsc_mat, SUBSET_NONZERO_PATTERN);
// }
