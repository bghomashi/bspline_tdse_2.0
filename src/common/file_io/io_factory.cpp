#include "common/file_io/io_factory.h"

using namespace io;

std::shared_ptr<IFactory> s_instance = nullptr;

bool IFactory::Startup() {
    return true;
}
void IFactory::Shutdown() {
}

ASCII IFactory::OpenASCII(const std::string& filename, char mode) {
    auto file = new IASCII();
    if (!file->Open(filename, mode))
        return nullptr;
    return ASCII(file);
}
void IFactory::CloseASCII(ASCII& o) {
    o = nullptr;
}


bool Factory::Startup() {
    return s_instance->Startup();
}
void Factory::Shutdown() {
    s_instance->Shutdown();
}

HDF5 Factory::OpenHDF5(const std::string& filename, char mode) {
    return s_instance->OpenHDF5(filename, mode);
}
void Factory::CloseHDF5(HDF5& o) {
    s_instance->CloseHDF5(o);
}

ASCII Factory::OpenASCII(const std::string& filename, char mode) {
    return s_instance->OpenASCII(filename, mode);
}
void Factory::CloseASCII(ASCII& o) {
    s_instance->CloseASCII(o);
}

void Factory::SetInstance(IFactory* factory) {
    s_instance = std::shared_ptr<IFactory>(factory);
}