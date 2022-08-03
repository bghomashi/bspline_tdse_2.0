#ifndef __PETSC_ASCI_H__
#define __PETSC_ASCI_H__

// class:       PetscASCII
// description: concrete implementation the "io_ascii" interface. Used the petsc/mpi
//              functions to ensure one the rank 0 process will interact with the file.

#include "common/file_io/io_ascii.h"
#include <string>

class PetscASCII : public io::IASCII {
public:
    PetscASCII();
    PetscASCII(const std::string& filename, char mode);
    ~PetscASCII();
 
    void Write(const std::string& text);
    std::string ReadLine();
    void Flush();
};

#endif