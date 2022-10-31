#include "petsc/file_io/petsc_binary.h"
#include "petsc/maths/petsc_common.h"
#include "common/utility/logger.h"
#include <cassert>
#include <petsc.h>

PetscBinary::PetscBinary() {
}
PetscBinary::PetscBinary(const std::string& filename, char mode) {
    PetscMPIInt rank;
    PetscErrorCode ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank);PETSCASSERT(ierr);
    if (rank == 0)
        if (!Open(filename, mode)) {
            Log::Critical("failed to open file: " + filename);
            PetscEnd();
        }
}
PetscBinary::~PetscBinary() {
    PetscMPIInt rank;
    PetscErrorCode ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank);PETSCASSERT(ierr);
    if (rank == 0)
        Close();
}
void PetscBinary::Write(const void* data, size_t size_in_bytes) {
    PetscMPIInt rank;
    PetscErrorCode ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank);PETSCASSERT(ierr);
    if (rank == 0) {
        IBinary::Write(data, size_in_bytes);
    }
}
void PetscBinary::Flush() {
    PetscMPIInt rank;
    PetscErrorCode ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank);PETSCASSERT(ierr);
    if (rank == 0) {
        IBinary::Flush();
    }
}
void PetscBinary::Read(void* data, size_t size_in_bytes) {
    PetscMPIInt rank;
    PetscErrorCode ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank);PETSCASSERT(ierr);
    std::string result;
    if (rank == 0)
        IBinary::Read(data, size_in_bytes);
}