#ifndef __CL_ENVIRONMENT_H__
#define __CL_ENVIRONMENT_H__

#define CL_USE_DEPRECATED_OPENCL_1_2_APIS
#include <CL/cl.h>
#include <string>

//************** Preprocessor directives ************* 
// 
// #define CL_AUTO_SYNC
// - defining to 1 will cause running kernels using the syntax 
//    "vector = kernel(args)"
//    to sync the vector before and after the run and all vector type arguments.
//    The can be a waste of bandwidth, pushing data that has not changed or 
//    reading data before it is used on the host side.
// #define CL_MANUAL_INIT
// - defining to 1 will cause the Environment (opencl context/queue/...) to
//   automatically initialize on first use.
//   Selects the default target device.

#ifndef CL_MANUAL_INIT 
#define CL_AUTO_INIT
#endif

namespace CL {
	typedef cl_double2 complex;
	class Environment {
		Environment();
	public:
		cl_device_id _device_id;
		cl_context _context;
		cl_command_queue _queue;
		size_t _max_wg, _local_mem, _global_mem;
		size_t _double_precision;

		cl_int _last_error;

		~Environment();
		static inline Environment& get() {
			static Environment instance;
			return instance;
		}
		inline cl_int getLastError() const {
			return _last_error;
		}
		inline size_t getMaxWorkGroups() const {
			return _max_wg;
		}
		inline size_t GetLocalMem() const {
			return _local_mem;
		}
		inline size_t GetGlobalMem() const {
			return _global_mem;
		}
		inline size_t GetDoublePrecisionSupport() const {
			return _double_precision;
		}

		// API - wrapper
		bool Initialize(unsigned target = 0);
		bool Destroy();

		// memory buffer api
		bool CreateBuffer(cl_mem& buf, cl_mem_flags flags, size_t size_in_bytes);
		bool DestroyBuffer(cl_mem buf);
		bool CreateSubBuffer(cl_mem buf, int from_byte, int size_in_bytes, cl_mem& sub_buffer);
		// kernel/program api
		bool EnqueueKernel(cl_kernel kernel, size_t dim, size_t* global_size, size_t* local_size);
		bool CreateProgramFromSource(cl_program& program, const std::string& source);
		bool DestroyProgram(cl_program program);
		bool DestroyKernel(cl_kernel kernel);
		bool GetKernelWorkGroupSize(cl_kernel& kernel, size_t& wg);
		bool CreateKernel(cl_kernel& kernel, cl_program program, const std::string& name);
		bool SetKernelArg(cl_kernel kernel, size_t n, cl_mem buf);
		bool SetKernelArg(cl_kernel kernel, size_t n, size_t size_in_bytes, const void* data);
		bool Read(cl_mem buf, size_t offset, size_t size_in_bytes, void* dest);
		bool Write(cl_mem buf, size_t offset, size_t size_in_bytes, const void* src);
		bool CopyBuffer(cl_mem src_buf, cl_mem dst_buf, size_t src_off, size_t dst_off, size_t size_in_bytes);
		void* Map(cl_mem buff, cl_map_flags flags, size_t offset, size_t bytes);
		void Unmap(cl_mem buff, void* host_ptr);

		// flush queue
		bool FlushQueue();
	};
}
#endif