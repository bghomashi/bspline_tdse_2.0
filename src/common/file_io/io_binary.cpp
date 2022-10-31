#include "common/file_io/io_binary.h"

using namespace io;

IBinary::IBinary() {

}
bool IBinary::Open(const std::string& filename, char mode) {
    switch (mode) {
        case 'w':
            _file = std::fstream(filename, std::ios::out | std::ios::binary);
            break;
        case 'r':
            _file = std::fstream(filename, std::ios::in | std::ios::binary);
            break;
        case 'a':
            _file = std::fstream(filename, std::ios::app | std::ios::binary);
            break;
    }
    return _file.is_open();
}
void IBinary::Close() {
    _file.close();
}

void IBinary::Write(const void* data, size_t size_in_bytes) {
    _file.write((const char*)data, size_in_bytes);
}
void IBinary::Flush() {
    _file.flush();
}
void IBinary::Read(void* data, size_t size_in_bytes) {
    _file.read((char*)data, size_in_bytes);
}

