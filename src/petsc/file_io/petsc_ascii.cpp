#include "petsc/file_io/petsc_ascii.h"
#include "petsc/maths/petsc_common.h"
#include "common/utility/logger.h"
#include <cassert>
#include <petsc.h>

PetscASCII::PetscASCII() {
}
PetscASCII::PetscASCII(const std::string& filename, char mode) {
    PetscMPIInt rank;
    PetscErrorCode ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank);PETSCASSERT(ierr);
    if (rank == 0)
        if (!Open(filename, mode)) {
            Log::Critical("failed to open file: " + filename);
            PetscEnd();
        }
}
PetscASCII::~PetscASCII() {
    PetscMPIInt rank;
    PetscErrorCode ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank);PETSCASSERT(ierr);
    if (rank == 0)
        Close();
}
void PetscASCII::Write(const std::string& text) {
    PetscMPIInt rank;
    PetscErrorCode ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank);PETSCASSERT(ierr);
    if (rank == 0) {
        IASCII::Write(text);
    }
}
void PetscASCII::Flush() {
    PetscMPIInt rank;
    PetscErrorCode ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank);PETSCASSERT(ierr);
    if (rank == 0) {
        IASCII::Flush();
    }
}
std::string PetscASCII::ReadLine() {
    PetscMPIInt rank;
    PetscErrorCode ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank);PETSCASSERT(ierr);
    std::string result;
    if (rank == 0)
        return IASCII::ReadLine();

    return result;
}