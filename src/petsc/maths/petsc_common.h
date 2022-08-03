#ifndef __PETSC_COMMON_H__
#define __PETSC_COMMON_H__

#include "petsc/maths/petsc_math_factory.h"
#include "petsc/maths/petsc_vector.h"
#include "petsc/maths/petsc_matrix.h"
#include "common/maths/vector.h"
#include <cassert>

#define PETSCASSERT(x) assert(x == 0)

PetscVector::Ptr_t PetscCast(maths::Vector o);
PetscMatrix::Ptr_t PetscCast(maths::Matrix o);

#endif