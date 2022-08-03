#include "petsc/maths/petsc_common.h"


std::shared_ptr<PetscVector> PetscCast(maths::Vector o) {
    return std::dynamic_pointer_cast<PetscVector>(o);
}
std::shared_ptr<PetscMatrix> PetscCast(maths::Matrix o) {
    return std::dynamic_pointer_cast<PetscMatrix>(o);
}