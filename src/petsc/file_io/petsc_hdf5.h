#ifndef __PETSC_HDF5_H__
#define __PETSC_HDF5_H__

#include "common/file_io/io_hdf5.h"
#include "common/maths/vector.h"
#include <petsc.h>

class PetscHDF5 : public io::IHDF5 {
    PetscViewer _viewer;
public:
    PetscHDF5();
    PetscHDF5(const std::string& filename, char mode);
    ~PetscHDF5();

    void PushGroup(const std::string& group_name);
    void PopGroup();
    bool HasGroup(const std::string& group_name) const;

    void ReadAttribute(const std::string& attr_name, maths::complex* value);
    void ReadAttribute(const std::string& attr_name, double* value);
    void ReadAttribute(const std::string& attr_name, int* value);
    void ReadAttribute(const maths::Vector object, const std::string& attr_name, maths::complex* value);
    void ReadAttribute(const maths::Vector object, const std::string& attr_name, double* value);
    void ReadAttribute(const maths::Vector object, const std::string& attr_name, int* value);
    void WriteAttribute(const maths::Vector object, const std::string& attr_name, const maths::complex value);
    void WriteAttribute(const maths::Vector object, const std::string& attr_name, const double value);
    void WriteAttribute(const maths::Vector object, const std::string& attr_name, const int value);
    void WriteAttribute(const std::string& attr_name, const maths::complex value);
    void WriteAttribute(const std::string& attr_name, const double value);
    void WriteAttribute(const std::string& attr_name, const int value);
    void ReadAttribute(const std::string& attr_name, const int value);
    bool HasAttribute(const std::string& attr_name) const;
    
    void WriteVector(const std::string& obj_name, const maths::Vector value);
    void ReadVector(const std::string& obj_name, maths::Vector value);
    bool HasVector(const std::string& obj_name) const;
};



#endif