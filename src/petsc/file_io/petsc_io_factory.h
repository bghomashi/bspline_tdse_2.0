#ifndef __IO_PETSC_FACTORY_H__
#define __IO_PETSC_FACTORY_H__

// class:       PetscFactory
// description: This class implements the io::Factory interface. It produces
//              io objects that use the PETSc API to only allow access to the 
//              resource on the rank 0 thread.

#include "common/file_io/io_factory.h"

class PetscIOFactory : public io::IFactory {
public:
    io::HDF5 OpenHDF5(const std::string& filename, char mode);
    void CloseHDF5(io::HDF5& o);

    io::ASCII OpenASCII(const std::string& filename, char mode);
    void CloseASCII(io::ASCII& o);
};

#endif