#ifndef __IO_BINARY_H__
#define __IO_BINARY_H__

// This class servers as a simple interface/wrapper for a binary 
// file. Implementing this interface allows the derived class
// to lock (mutex) the file resource but still be used by the rest 
// of the application without change, for example. 

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <memory>

namespace io {
    class IBinary {
    protected:
        std::fstream _file;
    public:
        IBinary();
        bool Open(const std::string& filename, char mode);
        void Close();

        virtual void Flush();
        virtual void Write(const void* data, size_t size_in_bytes);
        virtual void Read(void* data, size_t size_in_bytes);
    };
    
    typedef std::shared_ptr<IBinary> Binary;
}

#endif