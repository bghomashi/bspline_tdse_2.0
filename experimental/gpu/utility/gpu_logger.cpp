// #include "petsc/utility/petsc_logger.h"
// #include "petsc/maths/petsc_common.h"
// #include <petsc.h>

// void PetscLogger::Info(const std::string& text) {
//     PetscMPIInt rank;
//     PetscErrorCode ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank);PETSCASSERT(ierr);
//     if (rank == 0)
//         Logger::Info(text);
// }
// void PetscLogger::Warn(const std::string& text) {
//     PetscMPIInt rank;
//     PetscErrorCode ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank);PETSCASSERT(ierr);
//     if (rank == 0)
//         Logger::Warn(text);
// }
// void PetscLogger::Critical(const std::string& text) {
//     PetscMPIInt rank;
//     PetscErrorCode ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank);PETSCASSERT(ierr);
//     if (rank == 0)
//         Logger::Critical(text);
// }
// void PetscLogger::Debug(const std::string& text) {
//     PetscMPIInt rank;
//     PetscErrorCode ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank);PETSCASSERT(ierr);
//     if (rank == 0)
//         Logger::Debug(text);
// }
// void PetscLogger::SetLoggerFile(const std::string& log_file) {
//     PetscMPIInt rank;
//     PetscErrorCode ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank);PETSCASSERT(ierr);
//     if (rank == 0)
//         Logger::SetLoggerFile(log_file);
// }
// void PetscLogger::Flush() {
//     PetscMPIInt rank;
//     PetscErrorCode ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank);PETSCASSERT(ierr);
//     if (rank == 0)
//         Logger::Flush();
// }