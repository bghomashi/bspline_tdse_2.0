#include "cl_memory.hpp"


namespace CL {
	Buffer::Buffer() :
		e(Environment::get()), _dev_buf(0), _host_ptr(NULL), _bytes(0), _ref(new int(0)) {
	}
	Buffer::Buffer(Buffer& o) :
		e(Environment::get()), _dev_buf(o._dev_buf), _host_ptr(o._host_ptr), _bytes(o._bytes), _ref(o._ref) {
		(*_ref)++;
	}
	//Buffer::Buffer(Buffer&& o) :
	//	e(Environment::get()), _dev_buf(o._dev_buf), _host_ptr(o._host_ptr), _bytes(o._bytes), _ref(o._ref) {
	//	o._dev_buf = 0;
	//	o._host_ptr = 0;
	//	o._bytes = 0;
	//	(*_ref)++;
	//}
	Buffer::Buffer(size_t size_in_bytes, const void* init) :
		e(Environment::get()), _dev_buf(0), _host_ptr(NULL), _bytes(size_in_bytes), _ref(new int(0)) {
	
		cl_mem_flags flags = CL_MEM_READ_WRITE;
		if (!e.CreateBuffer(_dev_buf, flags, _bytes)) {
			_bytes = 0;
			return;
		}
		(*_ref)++;
		if (init != NULL)
			Write(init, 0, _bytes);
	}
	Buffer::~Buffer() {
		(*_ref)--;
		if (_dev_buf != 0 && *_ref == 0) {
			//Unmap();
			e.DestroyBuffer(_dev_buf);
		}
		_dev_buf = 0;
	}

	void Buffer::resize(int size_in_bytes) {
		(*this).~Buffer();
		new(this) Buffer(size_in_bytes);
	}
	void Buffer::Read(void* dest, size_t offset, size_t bytes) const {
		e.Read(_dev_buf, offset, bytes, dest);
	}
	void Buffer::Write(const void* src, size_t offset, size_t bytes) const {
		e.Write(_dev_buf, offset, bytes, src);
	}
	void Buffer::Copy(const Buffer& src, size_t src_offset, size_t dst_offset, size_t bytes) const {
		e.CopyBuffer(src._dev_buf, _dev_buf, src_offset, dst_offset, bytes);
	}
	void* Buffer::Map(cl_map_flags flags, size_t offset, size_t bytes) {
		if (_host_ptr != 0)
			Unmap();
		//_host_ptr = e.Map(_dev_buf, flags, offset, bytes == 0 ? _bytes : bytes);

		return _host_ptr;
	}
	void Buffer::Unmap() {
		if (_host_ptr != 0) {
			//e.Unmap(_dev_buf, _host_ptr);
			_host_ptr = 0;
		}
	}
	bool Buffer::CreateSubBuffer(int from_byte, int size_in_bytes, Buffer& buffer) {
		if (!e.CreateSubBuffer(_dev_buf, from_byte, size_in_bytes, buffer._dev_buf)) {
			return false;
		}

		_bytes = size_in_bytes;
		(*buffer._ref)++;

		return true;
	}
}