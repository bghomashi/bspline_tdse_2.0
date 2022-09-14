#include "gpu/maths/gpu_common.h"


std::shared_ptr<GPUVector> GPUCast(maths::Vector o) {
    return std::dynamic_pointer_cast<GPUVector>(o);
}
// std::shared_ptr<GPUMatrix> GPUCast(maths::Matrix o) {
//     return std::dynamic_pointer_cast<GPUMatrix>(o);
// }