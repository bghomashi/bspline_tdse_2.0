#include "petsc/file_io/petsc_io_factory.h"
#include "petsc/file_io/petsc_hdf5.h"
#include "petsc/file_io/petsc_ascii.h"

using namespace io;

HDF5 PetscIOFactory::OpenHDF5(const std::string& filename, char mode) {
    return HDF5(new PetscHDF5(filename, mode));
}
void PetscIOFactory::CloseHDF5(HDF5& file) {
    file = nullptr;
}


ASCII PetscIOFactory::OpenASCII(const std::string& filename, char mode) {
    return ASCII(new PetscASCII(filename, mode));
}
void PetscIOFactory::CloseASCII(ASCII& file) {
    file = nullptr;
}