#include "cl_kernel.hpp"


namespace CL {
	Program::Program() :
		e(Environment::get()), _program(0), _ref(new int(0)) {
	}
	Program::Program(const std::string& source) :
		e(Environment::get()), _program(0), _ref(new int(0)) {
		if (e.CreateProgramFromSource(_program, source))
			(*_ref)++;
	}
	Program::~Program() {
		(*_ref)--;
		if (_program != 0 && *_ref == 0)
			e.DestroyProgram(_program);
		_program = 0;
	}
	void Program::CompileFromSource(const std::string& source) {
		e.CreateProgramFromSource(_program, source);
	}

	Kernel Program::CreateKernel(const std::string& name, size_t dimensions) {
		Kernel k;
		if (e.CreateKernel(k._kernel, _program, name)) {
			(*k._ref)++;
			k._dimensions = dimensions;
			e.GetKernelWorkGroupSize(k._kernel, k._workgroup_size);
		}
		else
			k._kernel = 0;
		return k;
	}


	Kernel::Kernel() :
		e(Environment::get()), _kernel(0), _dimensions(0), _workgroup_size(0), _ref(new int(0)) {
	}
	Kernel::Kernel(const Kernel& kernel) noexcept :
		e(Environment::get()), _kernel(kernel._kernel), _dimensions(kernel._dimensions), _workgroup_size(kernel._workgroup_size), _ref(kernel._ref) {
		(*_ref)++;
	}
	Kernel::~Kernel() {
		(*_ref)--;
		if (_kernel != 0 && *_ref == 0)
			e.DestroyKernel(_kernel);
		_kernel = 0;
	}
	void Kernel::setArg(size_t n, const cl_mem& o) const {
		e.SetKernelArg(_kernel, n, o);
	}

	// recursively set the kernel arguments
	void Kernel::setArgs(size_t n) const {
		// done
	}

	// 
	void Kernel::Execute(size_t size) {
		e.EnqueueKernel(_kernel, 1, &size, NULL);
	}
	void Kernel::Run(size_t size) {
		if (_kernel == 0) return;
		size_t L[1] = { size };
		size_t s = std::min(_workgroup_size, e.getMaxWorkGroups());

		while (L[0] > s)
			L[0] = (L[0] + 1) / 2;

		size_t G[1] = { ((size + L[0] + 1) / L[0]) * L[0] };
	// std::cout << "global size = " << G[0] << std::endl;
	// std::cout << "local size = " << L[0] << std::endl;
		e.EnqueueKernel(_kernel, 1, G, L);
	}
	void Kernel::Run(size_t rows, size_t cols) {
		if (_kernel == 0) return;
		size_t L[2] = { rows, cols };
		size_t s = std::min(_workgroup_size, e.getMaxWorkGroups());

		while (L[0] > s)
			L[0] = (L[0] + 1) / 2;
		while (L[1] > s)
			L[1] = (L[1] + 1) / 2;


		while (L[0] * L[1] > s) {  // shrink the longer axis (square for locality)
			if (L[0] > 2 * L[1])
				L[0] = (L[0] + 1) / 2;
			else
				L[1] = (L[1] + 1) / 2;
		}

		size_t G[2] = { ((rows + L[0] - 1) / L[0]) * L[0], ((cols + L[1] - 1) / L[1]) * L[1] };

		e.EnqueueKernel(_kernel, 2, G, L);
	}

	std::string strip_parenthesis(const std::string& args) {
		std::string source = args;
		source = source.substr(1);
		source = source.substr(0, source.length() - 1);
		return source;
	}
	bool DumpSource(const std::string& filename, const std::string& source)
	{
		std::ofstream file(filename);
		if (file.is_open()) {
			file << source;
			return true;
		}
		return false;
	}
}