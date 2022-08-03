#ifndef __IO_ASCII_H__
#define __IO_ASCII_H__

// This class servers as a simple interface/wrapper for an ASCII 
// text file. Implementing this interface allows the derived class
// to lock (mutex) the file resource but still be used by the rest 
// of the application without change, for example. 

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <memory>

namespace io {
    class IASCII {
    protected:
        std::fstream _file;
    public:
        IASCII();
        bool Open(const std::string& filename, char mode);
        void Close();

        virtual void Flush();
        virtual void Write(const std::string& text);
        virtual std::string ReadLine();
    };
    
    typedef std::shared_ptr<IASCII> ASCII;
}

#endif