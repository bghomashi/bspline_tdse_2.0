#ifndef __IO_FACTORY_H__
#define __IO_FACTORY_H__

// class:       io::Factory
// description: This class is an interface to allow creating input/output
//              resources. Implentations of this interface might create 
//              io objects wrapped in a mutex.

#include "common/file_io/io_hdf5.h"
#include "common/file_io/io_ascii.h"
#include "common/file_io/io_binary.h"
#include <memory>

namespace io {
    class IFactory {
    public:
        virtual bool Startup();
        virtual void Shutdown();

        virtual HDF5 OpenHDF5(const std::string& filename, char mode) = 0;
        virtual void CloseHDF5(HDF5& o) = 0;

        virtual ASCII OpenASCII(const std::string& filename, char mode);
        virtual void CloseASCII(ASCII& o);

        virtual Binary OpenBinary(const std::string& filename, char mode);
        virtual void CloseBinary(Binary& o);
    };
    class Factory {
    public:
        static bool Startup();
        static void Shutdown();

        static HDF5 OpenHDF5(const std::string& filename, char mode);
        static void CloseHDF5(HDF5& o);

        static ASCII OpenASCII(const std::string& filename, char mode);
        static void CloseASCII(ASCII& o);

        static Binary OpenBinary(const std::string& filename, char mode);
        static void CloseBinary(Binary& o);

        static void SetInstance(IFactory* factory);
    };
}

#endif