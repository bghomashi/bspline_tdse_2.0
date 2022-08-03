#include "petsc/maths/petsc_math_factory.h"
#include "petsc/maths/petsc_vector.h"
#include "petsc/maths/petsc_matrix.h"
#include "petsc/maths/petsc_gmres.h"
#include "petsc/maths/petsc_eps.h"
#include "common/utility/logger.h"

using namespace maths;

bool PetscMathFactory::Startup(int argc, char **args) {
    PetscErrorCode ierr;
    ierr = PetscInitialize(&argc,&args,NULL,"Func Test\n"); 
    LOG_INFO("Initializing PETsc.");
    if (ierr) {
        LOG_CRITICAL("Failed!");
        return false;
    }
    ierr = SlepcInitialize(&argc,&args,NULL,NULL);if (ierr)
    if (ierr) {
        LOG_CRITICAL("Failed!");
        return false;
    }

    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &_rank); CHKERRQ(ierr);
    ierr = MPI_Comm_size(PETSC_COMM_WORLD,&_size); CHKERRQ(ierr);

    return true;
}
void PetscMathFactory::Shutdown() {
    PetscErrorCode ierr;
    ierr = SlepcFinalize();
    ierr = PetscFinalize();
}



Vector PetscMathFactory::CreateVector(int N) {
    return Vector(new PetscVector(N));
}
void PetscMathFactory::DestroyVector(Vector& v) {
    v = nullptr;                // If there are other references to v, the object is not destroyed
}

Matrix PetscMathFactory::CreateMatrix(int rows, int cols, int numBands) {
    return Matrix(new PetscMatrix(rows, cols, numBands));
}
void PetscMathFactory::DestroyMatrix(Matrix& m) {
    m = nullptr;                // If there are other references to m, the object is not destroyed
}
GMRESSolver PetscMathFactory::CreateGMRESSolver(int restart_iter, int max_iter) {
    return GMRESSolver(new PetscSolver(restart_iter, max_iter));
}
void PetscMathFactory::DestroyGMRESSolver(GMRESSolver& m) {
    m = nullptr;
}
maths::EigenSolver PetscMathFactory::CreateEigenSolver() {
    return maths::EigenSolver(new PetscEPS());
}
void PetscMathFactory::DestroyEigenSolver(maths::EigenSolver& m) {
    m = nullptr;
}