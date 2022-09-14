#ifndef __CL_MEMORY_H__
#define __CL_MEMORY_H__

#include <vector>
#include <memory>
#include <cstring>
#include "cl_environment.hpp"
#include "cl_kernel.hpp"

namespace CL {
	template <typename T>
	class Vector1D;
	template <typename T>
	class Vector2D;

	template<typename T>
	class Local {
	public:
		size_t n;
	};

	class Buffer {
	protected:
		std::shared_ptr<int> _ref;
	public:
		Environment& e;
		cl_mem _dev_buf;
		void* _host_ptr;
		size_t _bytes;
	public:
		Buffer();
		Buffer(Buffer& o);
		// Buffer(Buffer&& o);
		Buffer(size_t size_in_bytes, const void* init = NULL);
		~Buffer();

		void resize(int count);
		void Read(void* dest, size_t offset, size_t bytes) const;
		void Write(const void* src, size_t offset, size_t bytes) const;
		void Copy(const Buffer& src, size_t src_offset, size_t dst_offset, size_t bytes) const;
		void* Map(cl_map_flags flags = CL_MAP_READ + CL_MAP_WRITE, size_t offset = 0, size_t bytes = 0);
		void Unmap();
		bool CreateSubBuffer(int from_byte, int size_in_bytes, Buffer& buffer);

		inline operator cl_mem (void) const { return _dev_buf; }
		inline size_t get_byte_count(void) const { return _bytes; }
	};

	template <typename T>
	class Vector1D {
	public:
		int _count;
		std::vector<T> _host_data;
		Buffer _buffer;
	public:
		inline int size() const { return _count; }
		inline void shrink(int shrink_to) { _count = shrink_to; }
		inline operator cl_mem (void) const { return _buffer._dev_buf; }
		inline size_t get_byte_count(void) const { return _buffer._bytes; }
		void resize(int size_in_bytes);

		Vector1D();
		// Vector1D(Vector1D&& o);
		Vector1D(int count, const T* init = NULL);
		Vector1D(const std::initializer_list<T>& init);
		Vector1D(const std::vector<T>& init);

		void Read();
		void Write() const;
		void CopyHost(const Vector1D& src) ;
		void CopyGPU(const Vector1D& src) const;
		Vector1D CreateSubVector(int from, int size);
		bool CreateSubVector(int from, int count, Vector1D<T>& new_vector);

		inline std::vector<T>& host_vector() {
			return _host_data;
		}
		inline T* data() {
			return _host_data.data();
		}
		inline const T& operator() (size_t index) const {
			return _host_data[index];
		}
		inline T& operator() (size_t index) {
			return _host_data[index];
		}

		void operator=(Kernel& kernel) {
			if (kernel.getDimensions() != 1) {
				std::cout << "Kernel not properly vectorized.\n";
				return;
			}

	#ifdef CL_AUTO_SYNC
			Write();
	#endif
			kernel.setArg(0, _count);
			kernel.setArg(1, _buffer._dev_buf);
			kernel.Run(_count);

	#ifdef CL_AUTO_SYNC
			Read();
	#endif
		}
	};

	template <typename T>
	class Vector2D {
	public:
		int _rows, _cols;
		std::vector<T> _host_data;
		Buffer _buffer;
	public:
		inline int rows() const { return _rows; }
		inline int cols() const { return _cols; }
		inline int size() const { return _rows*_cols; }
		inline void shrink(int rows, int cols) { _rows = rows; _cols = cols; }

		Vector2D();
		Vector2D(int rows, int cols, const T* init = NULL);
		Vector2D(const std::initializer_list<std::initializer_list<T>>& init);
		Vector2D(int rows, int cols, const std::vector<T>& init);

		void Read();
		void Write() const;
		void CopyHost(const Vector2D& src) ;
		void CopyGPU(const Vector2D& src) const;

		inline const T& operator() (size_t row, size_t col) const {
			return _host_data[row + col*_rows];
		}
		inline T& operator() (size_t row, size_t col) {
			return _host_data[row + col*_rows];
		}

		void operator=(Kernel& kernel) {
			if (kernel.getDimensions() != 2) {
				std::cout << "Kernel not properly vectorized.\n";
				return;
			}

	#ifdef CL_AUTO_SYNC
			Write();
	#endif

			kernel.setArg(0, _rows);
			kernel.setArg(1, _cols);
			kernel.setArg(2, _buffer._dev_buf);
			kernel.Run(_rows, _cols);

	#ifdef CL_AUTO_SYNC
			Read();
	#endif
		}
	};

	
	template<typename T>
	Vector1D<T>::Vector1D()
		: _buffer(), _count(0) {}

	//template<typename T>
	//Vector1D<T>::Vector1D(Vector1D&& o) :
	//	e(Environment::get()), Buffer(o), _count(o._count), _host_data(o._host_data) {
	//	o._count = 0;
	//	o._host_data.clear();
	//}
	template<typename T>
	Vector1D<T>::Vector1D(int count, const T* init)
		: _buffer(count * sizeof(T), (const void*)init), _host_data(count), _count(count) {
		if (init != NULL)
			std::memcpy(&_host_data[0], init, count * sizeof(T));
	}
	template<typename T>
	Vector1D<T>::Vector1D(const std::initializer_list<T>& init)
		: _buffer(init.size() * sizeof(T), NULL), _host_data(init), _count(init.size()) {
			Write();
		}
	template<typename T>
	Vector1D<T>::Vector1D(const std::vector<T>& init)
		: _buffer(init.size() * sizeof(T), (void*)init.data()), _host_data(init), _count(init.size()) {}

	template<typename T>
	bool Vector1D<T>::CreateSubVector(int from, int count, Vector1D<T>& new_vector) {
		if (!_buffer.CreateSubBuffer(from*sizeof(T), count*sizeof(T), new_vector._buffer))
			return false;

		new_vector._count = count;
		// subvectors cannot sync data
		// new_vector._host_data.resize(0);

		return true;
	}
	template<typename T>
	void Vector1D<T>::resize(int count) {
		if (count <= _count) {
			shrink(count);
			return;
		}
		_buffer.resize(count * sizeof(T));
		_host_data.resize(count);
		_count = count;
	}
	template<typename T>
	void Vector1D<T>::Read() {
		_buffer.Read(&_host_data[0], 0, _count * sizeof(T));
	}
	template<typename T>
	void Vector1D<T>::Write() const {
		_buffer.Write(&_host_data[0], 0, _count * sizeof(T));
	}
	template<typename T>
	void Vector1D<T>::CopyHost(const Vector1D& src) {
		_host_data = src._host_data;
	}
	template<typename T>
	void Vector1D<T>::CopyGPU(const Vector1D& src) const {
		_buffer.Copy(src._buffer, 0, 0, _buffer.get_byte_count());
	}

	template<typename T>
	Vector2D<T>::Vector2D()
		: _buffer(), _rows(0), _cols(0) {}
	template<typename T>
	Vector2D<T>::Vector2D(int rows, int cols, const T* init)
		: _buffer(rows*cols*sizeof(T), (const void*)init), _host_data(rows*cols), _rows(rows), _cols(cols) {
		if (init != NULL)
			std::memcpy(&_host_data[0], init, rows*cols*sizeof(T));
	}
	template<typename T>
	Vector2D<T>::Vector2D(const std::initializer_list < std::initializer_list<T>>& init)
		: _rows(init.size()), _cols(init.begin().size()), _buffer(_rows*_cols*sizeof(T), NULL), _host_data(_rows*_cols, 0) {
		int i = 0, j = 0;
		for (const auto& a : init) {
			for (const auto& b : a) {
				_host_data[i + j*_rows] = b;
				++j;
			}
			j = 0;
			++i;
		}
		Write();
	}

	template<typename T>
	Vector2D<T>::Vector2D(int rows, int cols, const std::vector<T>& init)
		: _buffer(init.size()*sizeof(T), (void*)init.data()), _host_data(init), _rows(rows), _cols(cols) {}

	template<typename T>
	void Vector2D<T>::Read() {
		_buffer.Read(&_host_data[0], 0, _rows*_cols*sizeof(T));
	}
	template<typename T>
	void Vector2D<T>::Write() const {
		_buffer.Write(&_host_data[0], 0, _rows*_cols*sizeof(T));
	}
	template<typename T>
	void Vector2D<T>::CopyHost(const Vector2D& src) {
		_host_data = src._host_data;
	}
	template<typename T>
	void Vector2D<T>::CopyGPU(const Vector2D& src) const {
		_buffer.Copy(src._buffer, 0, 0, _buffer.get_byte_count());
	}
}
#endif