
#include "gpu/maths/gpu_vector.h"
#include "gpu/maths/gpu_common.h"
#include "gpu/kernels/cl_kernels.h"

using namespace maths;

GPUVector::GPUVector() {
    _len = 0;
}
GPUVector::GPUVector(int length) : _gpu_vec(length) {
    _len = length;
}
GPUVector::~GPUVector() {
}


complex GPUVector::Get(int index) const {
    return _gpu_vec(index);
}
void GPUVector::Set(int index, complex value) {
    _gpu_vec(index) = value;
}
void GPUVector::Scale(complex a) {
    // _gpu_vec = MathKernels::scale_vec(a);
}
void GPUVector::Duplicate(const Vector& o) {
    auto from = GPUCast(o);
    // destroy old vector if there was one
    _gpu_vec.resize(from->_len);
    _len = from->_len;
}
void GPUVector::Copy(const Vector& o) {
    auto from = GPUCast(o);
    _gpu_vec.CopyGPU(from->_gpu_vec);
}
void GPUVector::Zero() {
    // _gpu_vec = MathKernels::zero_vec()
}
double GPUVector::Max() const {
    // PetscReal max;
    // PetscErrorCode ierr = VecMax(_petsc_vec, NULL, &max); PETSCASSERT(ierr);
    return 0;
}
double GPUVector::Min() const {
    // PetscReal min;
    // PetscErrorCode ierr = VecMin(_petsc_vec, NULL, &min); PETSCASSERT(ierr);
    return 0;
}
void GPUVector::Abs() {
    // _gpu_vec = MathKernels::abs_vec();
}
void GPUVector::Concatenate(const std::vector<Vector>& vecs) {
    int nx = vecs.size();
    _len = 0;
    
    for (int i = 0; i < nx; i++)
        _len += GPUCast(vecs[i])->_len;

    _gpu_vec.resize(_len);

    int l = 0;
    for (int i = 0; i < nx; i++) {
        auto v = GPUCast(vecs[i]);
        _gpu_vec._buffer.Copy(
            v->_gpu_vec._buffer, 
            0, l*sizeof(complex), 
            v->_gpu_vec._buffer.get_byte_count());
        l += v->_len;
    }
}
void GPUVector::CopyTo(std::vector<complex>& values) {
    _gpu_vec.Read();
    _gpu_vec.CopyHost(values);
}

void GPUVector::Transform(Vector& out, std::function<std::vector<complex>(const std::vector<complex>&)> f) {
    _gpu_vec.Read();
    
    auto gpu_out = GPUCast(out);

    gpu_out->_gpu_vec.host_vector() = f(_gpu_vec.host_vector());

    gpu_out->_gpu_vec.Write();
}

Vector GPUVector::GetSubVector(int start, int end) {
    GPUVector* result = new GPUVector();

    result->_len = end - start;
    _gpu_vec.CreateSubVector(start, end - start, result->_gpu_vec);

    return Vector(result);
}
void GPUVector::RestoreSubVector(Vector sub) {
    auto gpu_sub = GPUCast(sub);
    gpu_sub->_gpu_vec.~Vector1D<complex>();
}


void GPUVector::Dot(const Vector b, complex& value) const {
    auto aa = this;
    auto bb = GPUCast(b);

    // Maths::VecDot(aa->_gpu_vec, bb->_gpu_vec, &value);
}
void GPUVector::AYPX(complex a, const Vector X) {
    auto y = this;
    auto x = GPUCast(X);

    // y->_gpu_vec = Maths::VecAYPX(a, x->_gpu_vec);
}
void GPUVector::AXPY(complex a, const Vector X) {
    auto y = this;
    auto x = GPUCast(X);

    // y->_gpu_vec = Maths::VecAXPY(a, x->_gpu_vec);
}



void GPUVector::AssembleBegin() {
}
void GPUVector::AssembleEnd() {
    _gpu_vec.Write();
}
