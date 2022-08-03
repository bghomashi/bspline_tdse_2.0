#ifndef __IO_HDF5_H__
#define __IO_HDF5_H__


// This class is an interface for an HDF5 file. 

#include <string>
#include <memory>

#include "common/maths/math_common.h"


namespace io {
    class IHDF5 {
    public:
        // ----- group management ----
        // "groups" in HDF5 files are like directories, so these functions allow 
        // the application to navigate the group structure. 
        virtual void PushGroup(const std::string& group_name) = 0;
        virtual void PopGroup() = 0;
        virtual bool HasGroup(const std::string& group_name) const = 0;

        // ----- read functions -----
        // allows the application to read an attribute or vector from the file
        // at the currect group
        virtual void ReadAttribute(const std::string& attr_name, maths::complex* value) = 0;
        virtual void ReadAttribute(const std::string& attr_name, double* value) = 0;
        virtual void ReadAttribute(const std::string& attr_name, int* value) = 0;
        virtual void ReadAttribute(const maths::Vector object, const std::string& attr_name, maths::complex* value) = 0;
        virtual void ReadAttribute(const maths::Vector object, const std::string& attr_name, double* value) = 0;
        virtual void ReadAttribute(const maths::Vector object, const std::string& attr_name, int* value) = 0;
        virtual void ReadVector(const std::string& obj_name, maths::Vector value) = 0;
        
        // ----- write functions ------
        // allows the application to write an attribute or vector from the file
        // at the currect group
        virtual void WriteAttribute(const maths::Vector object, const std::string& attr_name, const maths::complex value) = 0;
        virtual void WriteAttribute(const maths::Vector object, const std::string& attr_name, const double value) = 0;
        virtual void WriteAttribute(const maths::Vector object, const std::string& attr_name, const int value) = 0;
        virtual void WriteAttribute(const std::string& attr_name, const maths::complex value) = 0;
        virtual void WriteAttribute(const std::string& attr_name, const double value) = 0;
        virtual void WriteAttribute(const std::string& attr_name, const int value) = 0;
        virtual void WriteVector(const std::string& obj_name, const maths::Vector value) = 0;
        
        // ----- peek -----
        // check if an attribute or vector exists in the current group
        virtual bool HasAttribute(const std::string& attr_name) const = 0;
        virtual bool HasVector(const std::string& obj_name) const = 0;
    };

    typedef std::shared_ptr<IHDF5> HDF5;
}

#endif