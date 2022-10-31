#ifndef __PETSC_BINARY_H__
#define __PETSC_BINARY_H__

// class:       PetscBinary
// description: concrete implementation the "io_binary" interface. Used the petsc/mpi
//              functions to ensure one the rank 0 process will interact with the file.

#include "common/file_io/io_binary.h"
#include <string>

class PetscBinary : public io::IBinary {
public:
    PetscBinary();
    PetscBinary(const std::string& filename, char mode);
    ~PetscBinary();

    void Write(const void* data, size_t size_in_bytes);
    void Read(void* data, size_t size_in_bytes);
    void Flush();
};

#endif