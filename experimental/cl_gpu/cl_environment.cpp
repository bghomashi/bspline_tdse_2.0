#include "cl_environment.hpp"
#include <iostream>


namespace CL {
	Environment::Environment() :
		_device_id(NULL), _context(NULL), _queue(NULL), _last_error(0), _max_wg(0) {
#ifdef CL_AUTO_INIT
		Initialize();
#endif
	}
	Environment::~Environment() {
		Destroy();
	}

	bool Environment::Initialize(unsigned target) {
		if (_context != NULL) {
			std::cout << "Context already initialized! Cannot initialize twice.\n";
			return false;
		}
		const size_t MAX_PLAT = 8;
		cl_platform_id platforms[MAX_PLAT];
		cl_uint num_platforms = MAX_PLAT;
		_last_error = clGetPlatformIDs(MAX_PLAT, platforms, &num_platforms);
		if (_last_error != CL_SUCCESS) {
			std::cout << "Failed to get platform ids\n";
			return false;
		}
		cl_platform_id cpPlatform = platforms[0];
		cl_uint nDevices = 0;           // number of devices available
		cl_uint targetDevice = target;  // default device to compute on
		cl_uint nUnits;                 // number of compute units (SM's on NV GPU)

		int devtype = CL_DEVICE_TYPE_ALL;
		_last_error = clGetDeviceIDs(cpPlatform, devtype, 0, NULL, &nDevices);
		if (_last_error != CL_SUCCESS) {
			std::cout << "Failed to get device count\n";
			return false;
		}

		cl_device_id* devs = new cl_device_id[nDevices];
		std::cout << "number of devices: " << nDevices << std::endl;
		_last_error = clGetDeviceIDs(cpPlatform, devtype, nDevices, devs, NULL);
		if (_last_error != CL_SUCCESS) {
			std::cout << "Failed to get device ids\n";
			return false;
		}
		_device_id = devs[targetDevice];
		delete[] devs;

		_last_error = clGetDeviceInfo(_device_id, CL_DEVICE_MAX_COMPUTE_UNITS, sizeof(nUnits), &nUnits, NULL);
		if (_last_error != CL_SUCCESS) {
			std::cout << "Failed to get device info\n";
			return false;
		}

		// this can be used to allow interop between OpenGL buffers
		cl_context_properties* props = NULL;
		_context = clCreateContext(NULL, 1, &_device_id, NULL, NULL, &_last_error);
		if (_last_error != CL_SUCCESS) {
			std::cout << "Failed to get create context\n";
			return false;
		}
		_queue = clCreateCommandQueue(_context, _device_id, 0, &_last_error);
		if (_last_error != CL_SUCCESS) {
			std::cout << "Failed to get create command queue\n";
			return false;
		}

		_last_error = clGetDeviceInfo(_device_id, CL_DEVICE_MAX_WORK_GROUP_SIZE, sizeof(_max_wg), &_max_wg, 0);
		_last_error = clGetDeviceInfo(_device_id, CL_DEVICE_LOCAL_MEM_SIZE, sizeof(_local_mem), &_local_mem, 0);
		_last_error = clGetDeviceInfo(_device_id, CL_DEVICE_DOUBLE_FP_CONFIG, sizeof(_double_precision), &_double_precision, 0);
		_last_error = clGetDeviceInfo(_device_id, CL_DEVICE_GLOBAL_MEM_SIZE, sizeof(_global_mem), &_global_mem, 0);

		std::cout << "OpenCL initialized with " << nUnits << " compute units." << std::endl;

		return true;
	}


	bool Environment::CreateBuffer(cl_mem& buf, cl_mem_flags flags, size_t size_in_bytes) {
		buf = clCreateBuffer(_context, flags, size_in_bytes, NULL, &_last_error);
		switch (_last_error) {
		case CL_INVALID_CONTEXT:
			std::cout << "CreateBuffer: Not a valid context." << std::endl;
			return false;
		case CL_INVALID_VALUE:
			std::cout << "CreateBuffer: Invalid flags." << std::endl;
			return false;
		case CL_INVALID_BUFFER_SIZE:
			std::cout << "CreateBuffer: Buffer size must be greater than 0 and less	than CL_DEVICE_MAX_MEM_ALLOC_SIZE ." << std::endl;
			return false;
		case CL_INVALID_HOST_PTR:					// this is current assumed to be NULL. Not sure what it is for.
			//if (host_ptr == NULL)
			std::cout << "CreateBuffer: host_ptr is NULL and (CL_MEM_USE_HOST_PTR || CL_MEM_COPY_HOST_PTR)." << std::endl;
			//else
			//	std::cout << "CreateBuffer: host_ptr is not NULL and not (CL_MEM_COPY_HOST_PTR || CL_MEM_USE_HOST_PTR) OR  host." << std::endl;
			return false;
		case CL_MEM_OBJECT_ALLOCATION_FAILURE:
			std::cout << "CreateBuffer: Failer to allocate buffer memory (GPU side)." << std::endl;
			return false;
		case CL_OUT_OF_HOST_MEMORY:
			std::cout << "CreateBuffer: Failure to allocate memory on the host (CPU side)." << std::endl;
			return false;
		}
		return true;
	}
	bool Environment::DestroyBuffer(cl_mem buf) {
		_last_error = clReleaseMemObject(buf);
		return true;
	}

	bool Environment::EnqueueKernel(cl_kernel kernel, size_t dim, size_t* global_size, size_t* local_size) {
		_last_error = clEnqueueNDRangeKernel(_queue, kernel, dim, NULL, global_size, local_size, 0, NULL, NULL);
		if (_last_error != CL_SUCCESS) {
			std::cout << "EnqueueKernel: failed to execute kernel." << std::endl;
			return false;
		}
		return true;
	}

	bool Environment::Read(cl_mem buf, size_t offset, size_t size_in_bytes, void* dest) {
		_last_error = clEnqueueReadBuffer(_queue, buf, CL_TRUE, offset, size_in_bytes, dest, 0, NULL, NULL);
		if (_last_error != CL_SUCCESS) {
			std::cout << "Read: failed to read memory buffer." << std::endl;
			return false;
		}
		return true;
	}

	bool Environment::Write(cl_mem buf, size_t offset, size_t size_in_bytes, const void* src) {
		_last_error = clEnqueueWriteBuffer(_queue, buf, CL_TRUE, offset, size_in_bytes, src, 0, NULL, NULL);
		if (_last_error != CL_SUCCESS) {
			std::cout << "Write: failed to write memory buffer." << std::endl;
			return false;
		}

		return true;
	}
	bool Environment::CopyBuffer(cl_mem src_buf, cl_mem dst_buf, size_t src_off, size_t dst_off, size_t size_in_bytes) {
		_last_error = clEnqueueCopyBuffer(_queue, src_buf, dst_buf, src_off, dst_off, size_in_bytes, 0, NULL, NULL);
		if (_last_error != CL_SUCCESS) {
			std::cout << "CopyBuffer: failed to copy memory buffer." << std::endl;
			return false;
		}

		return true;
	}
	
	bool Environment::CreateSubBuffer(cl_mem buf, int from_byte, int size_in_bytes, cl_mem& sub_buffer) {
    	cl_buffer_region region;
		region.origin = from_byte;
    	region.size = size_in_bytes;
		sub_buffer = clCreateSubBuffer(buf, CL_MEM_READ_ONLY,
                               CL_BUFFER_CREATE_TYPE_REGION, &region, &_last_error);
		if (_last_error != CL_SUCCESS) {
			std::cout << "CreateSubBuffer: failed to create sub-buffer." << std::endl;
			return false;
		}

		return true;
	}
	void* Environment::Map(cl_mem buff, cl_map_flags flags, size_t offset, size_t bytes) {
		void* host_ptr = clEnqueueMapBuffer(_queue, buff,		// queue, buffer
			CL_TRUE, flags,									    // block?, flags
			offset, bytes,			                            // offset, bytes
			0, 0, 0, &_last_error);								// events
		if (_last_error != CL_SUCCESS) {
			std::cout << "Map: Failed to map memory buffer." << std::endl;
		}
		return host_ptr;
	}
	void Environment::Unmap(cl_mem buff, void* host_ptr) {

		_last_error = clEnqueueUnmapMemObject(_queue, buff,
			host_ptr,
			0, 0, 0);
		if (_last_error != CL_SUCCESS) {
			std::cout << "Unmap: Failed to unmap memory buffer." << std::endl;
		}
	}

	bool Environment::CreateProgramFromSource(cl_program& program, const std::string& source) {
		const char* source_str = source.c_str();
		size_t source_size = source.length();
		program = clCreateProgramWithSource(_context, 1, (const char**)&source_str, (const size_t*)&source_size, &_last_error);
		if (_last_error != CL_SUCCESS) {
			std::cout << "CreateProgramFromSource: failed to create program from source." << std::endl;
			return false;
		}
		// Build the program
		_last_error = clBuildProgram(program, 1, &_device_id, NULL, NULL, NULL);
		if (_last_error != CL_SUCCESS) {
			std::cout << "CreateProgramFromSource: failed to build program." << std::endl;

			size_t logsize = 0;
			_last_error = clGetProgramBuildInfo(program, _device_id, CL_PROGRAM_BUILD_LOG, 0, NULL, &logsize);
			char* build_log = new char[logsize + 1];
			_last_error = clGetProgramBuildInfo(program, _device_id, CL_PROGRAM_BUILD_LOG, logsize, build_log, NULL);
			build_log[logsize] = 0;
			std::cout << "Build log: " << build_log << std::endl;
			delete[] build_log;

			return false;
		}
		return true;
	}

	bool Environment::CreateKernel(cl_kernel& kernel, cl_program program, const std::string& name) {
		kernel = clCreateKernel(program, name.c_str(), &_last_error);
		if (_last_error != CL_SUCCESS) {
			std::cout << "CreateKernel: failed to create kernel." << std::endl;
			return false;
		}
		return true;
	}

	bool Environment::GetKernelWorkGroupSize(cl_kernel& kernel, size_t& wg) {
		_last_error = clGetKernelWorkGroupInfo(kernel, _device_id, CL_KERNEL_WORK_GROUP_SIZE, sizeof(size_t), &wg, 0);
		if (_last_error != CL_SUCCESS) {
			std::cout << "GetKernelWorkGroupSize: failed to get kernel workgroup size." << std::endl;
			return false;
		}
		return true;
	}
	bool Environment::SetKernelArg(cl_kernel kernel, size_t n, cl_mem buf) {
		_last_error = clSetKernelArg(kernel, n, sizeof(cl_mem), (void*)&buf);
		if (_last_error != CL_SUCCESS) {
			std::cout << "SetKernelArg: failed to set kernel argument." << std::endl;
			return false;
		}
		return true;
	}
	bool Environment::SetKernelArg(cl_kernel kernel, size_t n, size_t size_in_bytes, const void* data) {
		_last_error = clSetKernelArg(kernel, n, size_in_bytes, data);
		if (_last_error != CL_SUCCESS) {
			std::cout << "SetKernelArg: failed to set kernel argument." << std::endl;
			return false;
		}
		return true;
	}

	bool Environment::DestroyProgram(cl_program program) {
		_last_error = clReleaseProgram(program);
		return true;
	}
	bool Environment::DestroyKernel(cl_kernel kernel) {
		_last_error = clReleaseKernel(kernel);
		return true;
	}
	bool Environment::Destroy() {
		if (_queue != 0)
			_last_error = clReleaseCommandQueue(_queue);
		if (_context != 0)
			_last_error = clReleaseContext(_context);

		_queue = 0;
		_context = 0;

		return true;
	}
	bool Environment::FlushQueue() {
		_last_error = clFlush(_queue);
		_last_error = clFinish(_queue);
		return true;
	}
}