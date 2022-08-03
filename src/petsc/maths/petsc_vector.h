#ifndef __PETSC_VECTOR_H__
#define __PETSC_VECTOR_H__

#include "common/maths/math_common.h"
#include "common/maths/vector.h"
#include <cassert>
#include <petsc.h>
#include <slepc.h>

#define PETSCASSERT(x) assert(x == 0)

class PetscEPS;
class PetscGMRES;
class PetscEngine;

class PetscVector : public maths::IVector {
    friend PetscEngine;
    friend PetscEPS;
    friend PetscGMRES;

    friend void MatMult(const maths::Matrix& M, const maths::Vector& in, maths::Vector& out);
    friend void MatAXPY(maths::Matrix& X, maths::complex a, const maths::Matrix& Y);

    bool _dirty;
    std::vector<maths::complex> _local_copy;
public:
    typedef std::shared_ptr<PetscVector> Ptr_t;

    VecScatter _petsc_ctx;
    Vec _petsc_sca_vec;
    Vec _petsc_vec;
    IS _petsc_is;

    PetscVector();
    PetscVector(int length);
    ~PetscVector();

    // filling in the maths::vector API
    maths::complex Get(int index) const;
    void Set(int index, maths::complex value);
    void Scale(maths::complex a);
    void Duplicate(const maths::Vector& o);
    void Copy(const maths::Vector& o);
    void Zero();
    void Concatenate(const std::vector<maths::Vector>& vecs);
    void CopyTo(std::vector<maths::complex>& values); 
    void Transform(maths::Vector& out, std::function<std::vector<maths::complex>(const std::vector<maths::complex>&)> f);
    maths::Vector GetSubVector(int start, int end);
    void RestoreSubVector(maths::Vector sub);

    void AYPX(maths::complex a, const maths::Vector X);
    void AXPY(maths::complex a, const maths::Vector X);
    void Dot(const maths::Vector b, maths::complex& value) const;


    void AssembleBegin();
    void AssembleEnd();

    void CreateScatter();
    void Scatter();
    void ScatterGetArray(maths::complex** ptr);
    void ScatterRestoreArray(maths::complex** ptr);
};

#endif