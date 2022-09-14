#ifndef __CL_KERNEL_H__
#define __CL_KERNEL_H__

#include <iostream>
#include <memory>
#include <fstream>
#include <string>
#include "cl_environment.hpp"


namespace CL {
	class Kernel;
	class Buffer;
	template<typename T>
	class Vector1D;
	template<typename T>
	class Vector2D;
	template<typename T>
	class Local;

	class Program {
		std::shared_ptr<int> _ref;
	public:
		Environment& e;
		cl_program _program;
	public:
		Program();
		Program(const std::string& source);
		~Program();
		Program& operator=(const Program& o) {
			_program = o._program;
			_ref = o._ref;
			(*_ref)++;
			return *this;
		}

		void CompileFromSource(const std::string& source);
		Kernel CreateKernel(const std::string& name, size_t dimensions = 0);
	};

	class Kernel {
		friend Program;
		std::shared_ptr<int> _ref;
	public:
		Environment& e;

		cl_kernel _kernel;
		size_t _workgroup_size;
		size_t _dimensions;
	public:
		Kernel();
		Kernel(const Kernel& kernel) noexcept;
		~Kernel();

		inline size_t getDimensions() const { return _dimensions; }
		inline void setDimensions(size_t dimensions) { _dimensions = dimensions; }
		inline size_t getWorkGroupSize() const {
			return  std::min(_workgroup_size, e.getMaxWorkGroups());
		}
		// 
		void Run(size_t size);
		void Run(size_t rows, size_t cols);
		void Execute(size_t size);

		inline Kernel& operator=(const Kernel& o) {
			_kernel = o._kernel;
			_workgroup_size = o._workgroup_size;
			_dimensions = o._dimensions;
			_ref = o._ref;
			(*_ref)++;

			return *this;
		}

		// allows syntax kernel(args)
		template <typename ... Ts>
		Kernel& operator() (Ts... args);

		// set kernel arguments
		// overload for cl_mem object directly
		void setArg(size_t n, const cl_mem& o) const;
		// overload for local buffer
		template<typename T>
		void setArg(size_t n, const Local<T>& o) const;
		// overload for Vector1D object directly
		template<typename T>
		void setArg(size_t n, const Vector1D<T>& o) const;
		// overload for Vector2D object directly
		template<typename T>
		void setArg(size_t n, const Vector2D<T>& o) const;
		// overload for Buffer object
		template <typename T> typename std::enable_if<std::is_base_of<Buffer, T>::value>::type
			setArg(size_t n, const T& o) const;
		// catch all overload for everything else
		template <typename T> typename std::enable_if<!std::is_base_of<Buffer, T>::value>::type
			setArg(size_t n, const T& o) const;

	private:
		// recursive functions to allow () operator to function
		void setArgs(size_t n) const;
		template <typename T, typename ... Ts>
		void setArgs(size_t n, const T& t, Ts... args) const;
	};

	template <typename ... Ts>
	Kernel& Kernel::operator() (Ts... args) {
		size_t start_n = _dimensions == 0 ? 0 : _dimensions + 1;
		setArgs(start_n, args...);
		return *this;
	}
	template<typename T>
	void Kernel::setArg(size_t n, const Local<T>& o) const {
		e.SetKernelArg(_kernel, n, sizeof(T)*o.n, NULL);
	}
	template<typename T>
	void Kernel::setArg(size_t n, const Vector1D<T>& o) const {
#ifdef CL_AUTO_SYNC
		o.Write();
#endif
		e.SetKernelArg(_kernel, n, o._buffer._dev_buf);
	}
	template<typename T>
	void Kernel::setArg(size_t n, const Vector2D<T>& o) const {
#ifdef CL_AUTO_SYNC
		o.Write();
#endif
		e.SetKernelArg(_kernel, n, o._buffer._dev_buf);
	}
	template <typename T> typename std::enable_if<std::is_base_of<Buffer, T>::value>::type
	Kernel::setArg(size_t n, const T& o) const {
		e.SetKernelArg(_kernel, n, o._dev_buf);
	}
	template <typename T> typename std::enable_if<!std::is_base_of<Buffer, T>::value>::type
	Kernel::setArg(size_t n, const T& o) const {
		e.SetKernelArg(_kernel, n, sizeof(T), (const void*)&o);
	}
	template <typename T, typename ... Ts>
	void Kernel::setArgs(size_t n, const T& t, Ts... args) const {
		setArg(n, t);
		setArgs(n + 1, args...);
	}
	std::string strip_parenthesis(const std::string& args);
}
	// define __global as nothing to get rid of intelisense errors
	#define __global

	/****** Macro to assist in source for 1D kernel ************
	* global variables:
	*     i		- index of this kernel instance
	*     size	- total number of instances
	*     result - variable to store instance result. Also contains incoming value from Vector1D
	*
	* usage example:
	* std::string source = Kernel1D(
	*     int, vector_add, (__global const int* A, __global const int* B),
	*         result = A[i] + B[i];
	*     );
	**************************************************************/

	#define Kernel1D(ret_type, name, args, op) \
		"\n\n__kernel void " #name "(int size, __global " #ret_type " *ret_array, " + CL::strip_parenthesis(#args) + ") {\n" \
		"    int i = get_global_id(0);\n" \
		"    if (i < size) {\n" \
		"        " #ret_type " result = ret_array[i];\n" \
		"        // Do the operation \n" \
		"        " #op ";\n" \
		"        ret_array[i] = result;\n" \
		"    }\n" \
		"}\n"

	/****** Macro to assist in source for 1D kernel ************
	* global variables:
	*     idx	- index of this kernel instance
	*     r     - row of this kernel instance
	*     c     - row of this kernel instance
	*     rows	- total number of rows
	*     cols	- total number of cols
	*     result - variable to store instance result. Also contains incoming value from Vector2D
	*
	* usage example:
	* std::string source = Kernel2D(
	*     int, matrix_add, (__global const int* A, __global const int* B, __global const int* C),
	*         result = A[r] + B[c] + C[idx];
	*     );
	**************************************************************/
	#define Kernel2D(ret_type, name, args, op) \
		"\n\n__kernel void " #name "(int rows, int cols, __global " #ret_type " *ret_array, " + CL::strip_parenthesis(#args) + ") {\n" \
		"    int r = get_global_id(0);\n" \
		"    int c = get_global_id(1);\n" \
		"    if (r < rows && c < cols) {\n" \
		"        const int idx = r + c*rows;\n" \
		"        "#ret_type" result = ret_array[idx];\n" \
		"        " #op ";\n" \
		"        ret_array[idx] = result;\n" \
		"    }\n" \
		"}\n"

	#define StaticCode(src) "\n\n" #src "\n\n"
#endif