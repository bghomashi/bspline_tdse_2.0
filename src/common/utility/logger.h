#ifndef __LOGGER_H__
#define __LOGGER_H__

// class:       Logger
// description: The Logger class/interface implements a basic logger with several 
//              "log levels". The purpose is to allow the user to set the desired
//              level at compile time and log calls above this levelwill be 
//              reported.
//
//              ** Do not use the class directly. Instead use the static API 
//              defined at the start/end of the file **



#include <string>
#include <memory>

#define LOG_LEVEL_NONE        0
#define LOG_LEVEL_INFO        1
#define LOG_LEVEL_WARN        2
#define LOG_LEVEL_CRITICAL    3
#define LOG_LEVEL_DEBUG       4

#ifndef LOG_LEVEL
#define LOG_LEVEL LOG_LEVEL_DEBUG
#endif

class Logger {
public:
    typedef std::shared_ptr<Logger> Ptr_t;
    
    // write text to the various logger levels
    virtual void Info(const std::string& text);
    virtual void Warn(const std::string& text);
    virtual void Critical(const std::string& text);
    virtual void Debug(const std::string& text);

    // set the logger to output to a file
    virtual void SetLoggerFile(const std::string& log_file);
    // flush the log buffer
    virtual void Flush();
};


// static API for the logger. If a derived class is implemented set the
// logger to an instance of the derived class:
//
// Log::SetLogger(new DerivedLogger());

namespace Log {
    void Info(const std::string& text);
    void Warn(const std::string& text);
    void Critical(const std::string& text);
    void Debug(const std::string& text);
    void SetLoggerFile(const std::string& log_file);
    void SetLogger(Logger* logger);
    void Flush();
}

#if LOG_LEVEL >= LOG_LEVEL_INFO
#define LOG_INFO(x) Log::Info(x)
#else
#define LOG_INFO(x)
#endif

#if LOG_LEVEL >= LOG_LEVEL_WARN
#define LOG_WARN(x) Log::Warn(x)
#else
#define LOG_WARN(x)
#endif

#if LOG_LEVEL >= LOG_LEVEL_CRITICAL
#define LOG_CRITICAL(x) Log::Critical(x)
#else
#define LOG_CRITICAL(x)
#endif

#if LOG_LEVEL >= LOG_LEVEL_DEBUG
#define LOG_DEBUG(x) Log::Debug(x)
#else
#define LOG_DEBUG(x)
#endif

#define LOG_FLUSH Log::Flush()

#endif