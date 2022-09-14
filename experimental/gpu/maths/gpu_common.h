#ifndef __GPU_COMMON_H__
#define __GPU_COMMON_H__

// #include "gpu/maths/gpu_math_factory.h"
#include "gpu/maths/gpu_vector.h"
// #include "gpu/maths/gpu_matrix.h"
#include "common/maths/vector.h"
#include <cassert>

GPUVector::Ptr_t GPUCast(maths::Vector o);
// GPUMatrix::Ptr_t GPUCast(maths::Matrix o);

#endif