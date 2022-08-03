#ifndef __PETSC_LOGGER_H__
#define __PETSC_LOGGER_H__

// class:       PetscLogger
// description: Implements the logger interface. The purpose of this class
//              is to restrict access to the log resource to only the rank 0
//              process.


#include "common/utility/logger.h"
#include <string>

class PetscLogger : public Logger {
public:
    void Info(const std::string& text);
    void Warn(const std::string& text);
    void Critical(const std::string& text);
    void Debug(const std::string& text);
    void SetLoggerFile(const std::string& log_file);
    void Flush();
};



#endif