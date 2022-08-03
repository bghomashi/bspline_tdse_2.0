#include <iostream>
#include <fstream>
#include "common/utility/logger.h"

static Logger::Ptr_t s_logger(new Logger());
static std::ofstream s_logger_file;

void Logger::Info(const std::string& text) {
    if (s_logger_file && s_logger_file.good())
        s_logger_file << "INFO: " << text << std::endl;
    else
        std::cout << "INFO: " << text << std::endl;
}
void Logger::Warn(const std::string& text) {
    if (s_logger_file && s_logger_file.good())
        s_logger_file << "WARN: " << text << std::endl;
    else
        std::cout << "WARN: " << text << std::endl;
}
void Logger::Critical(const std::string& text) {
    if (s_logger_file && s_logger_file.good())
        s_logger_file << "CRITICAL: " << text << std::endl;
    else
        std::cout << "CRITICAL: " << text << std::endl;
}
void Logger::Debug(const std::string& text) {
    if (s_logger_file && s_logger_file.good())
        s_logger_file << "DEBUG: " <<  text << std::endl;
    else
        std::cout << "DEBUG: " << text << std::endl;
}
void Logger::SetLoggerFile(const std::string& log_file) {
    s_logger_file = std::ofstream(log_file, std::ostream::app);
}
void Logger::Flush() {
    if (s_logger_file && s_logger_file.good())
        s_logger_file << std::flush;
    else
        std::cout << std::flush;
}

// these are the static defined functions
namespace Log {
    void Info(const std::string& text) {
        s_logger->Info(text);
    }
    void Warn(const std::string& text) {
        s_logger->Warn(text);
    }
    void Critical(const std::string& text) {
        s_logger->Critical(text);
    }
    void Debug(const std::string& text) {
        s_logger->Debug(text);
    }
    void SetLoggerFile(const std::string& log_file) {
        s_logger->SetLoggerFile(log_file);
    }
    void SetLogger(Logger* logger) {
        s_logger.reset(logger);
    }
    void Flush() {
        s_logger->Flush();
    }
}