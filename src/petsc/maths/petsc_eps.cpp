
#include "petsc/maths/petsc_eps.h"
#include "petsc/maths/petsc_matrix.h"
#include "petsc/maths/petsc_common.h"
#include <petsc.h>
#include <slepc.h>

using namespace maths;

void PetscEPS::Solve(const Matrix A, const Matrix S, int numVectors, double tol) {
    PetscErrorCode ierr;
    
    complex norm;
    Vec tempDot;
    EPS petsc_eps;
    int nconv;
    int rank;

    // cast matrices to have access to petsc variables
    PetscMatrix::Ptr_t petscS = PetscCast(S), 
                       petscA = PetscCast(A);

    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); PETSCASSERT(ierr);
    ierr = EPSCreate( PETSC_COMM_WORLD, &petsc_eps ); PETSCASSERT(ierr);
    ierr = EPSSetOperators( petsc_eps, petscA->_petsc_mat, petscS->_petsc_mat ); PETSCASSERT(ierr);
    ierr = EPSSetTolerances(petsc_eps, tol, PETSC_DEFAULT); PETSCASSERT(ierr);
    switch (_problemType)
    {
    case EigenProblemType::GHEP:
        ierr = EPSSetProblemType(petsc_eps, EPS_GHEP); PETSCASSERT(ierr);
        break;
    case EigenProblemType::GNHEP:
        ierr = EPSSetProblemType(petsc_eps, EPS_GNHEP); PETSCASSERT(ierr);
        break;
    case EigenProblemType::PGNHEP:
        ierr = EPSSetProblemType(petsc_eps, EPS_PGNHEP); PETSCASSERT(ierr);
        break;
    }
    
    ierr = EPSSetDimensions(petsc_eps, numVectors, PETSC_DEFAULT, std::max(petscA->Rows(),100)); PETSCASSERT(ierr);
    ierr = EPSSetWhichEigenpairs(petsc_eps, EPS_SMALLEST_REAL); PETSCASSERT(ierr);

    ierr = EPSSetFromOptions(petsc_eps); PETSCASSERT(ierr);

    ierr = EPSSolve(petsc_eps); PETSCASSERT(ierr);
    ierr = EPSGetConverged(petsc_eps, &nconv); PETSCASSERT(ierr);

    _values.clear();
    _vectors.clear();
    _values.resize(std::min(nconv-1, numVectors));
    _vectors.reserve(std::min(nconv-1, numVectors));
    
    // make a vector just to hold intermediate dot-product later
    MatCreateVecs(petscS->_petsc_mat, &tempDot, NULL);

    for (int j = 0; j < _values.size(); j++) {
        PetscVector* temp = new PetscVector();

        MatCreateVecs(petscS->_petsc_mat, &temp->_petsc_vec, NULL);
        temp->_len = petscS->Rows();
        
        // get the eigen-pair
        EPSGetEigenpair(petsc_eps, j, &_values[j], NULL, temp->_petsc_vec, NULL);
        
        // ensure eigenvector is real
        complex sign = 0.0;
        const complex* xx;
        if (rank == 0) {
            VecGetArrayRead(temp->_petsc_vec, &xx);
            sign = *xx/std::abs(*xx);                   // grab overall phase
            VecRestoreArrayRead(temp->_petsc_vec, &xx);
        }
        MPI_Bcast(&sign, 1, MPIU_SCALAR, 0, PETSC_COMM_WORLD);
        VecScale(temp->_petsc_vec, 1./sign);            // remove phase
        VecRealPart(temp->_petsc_vec);

        // check sign convention
        if (rank == 0) {
            VecGetArrayRead(temp->_petsc_vec, &xx);
            sign = *xx/std::abs(*xx);                   // ensure the first coefficient is positive
            VecRestoreArrayRead(temp->_petsc_vec, &xx);
        }
        MPI_Bcast(&sign, 1, MPIU_SCALAR, 0, PETSC_COMM_WORLD);
        VecScale(temp->_petsc_vec, sign);            // remove phase


        // normalize
        MatMult(petscS->_petsc_mat, temp->_petsc_vec, tempDot);
        VecDot(temp->_petsc_vec, tempDot, &norm);
        VecScale(temp->_petsc_vec, 1./sqrt(norm));

        // MatMult(petscS->_petsc_mat, temp->_petsc_vec, tempDot);
        // VecDot(temp->_petsc_vec, tempDot, &norm);
    
        _vectors.push_back(Vector(temp));
        
        //if (_rank == 0)
            // PetscPrintf(PETSC_COMM_SELF, "energy %d: (%f, %f)\n", j, std::real(values[j]), std::imag(values[j]));
            // PetscPrintf(PETSC_COMM_SELF, "energy %d: (%f, %f) and norm: %f\n", j, std::real(values[j]), std::imag(values[j]), std::real(norm));
    }

    EPSDestroy(&petsc_eps);
    VecDestroy(&tempDot);
}