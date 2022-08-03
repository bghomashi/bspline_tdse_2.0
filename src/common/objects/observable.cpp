#include "common/objects/observable.h"

using namespace tdse;

Observable::Observable() : _computePeriod(1) {
}
void Observable::DoObservable(int it) {
    if (it % _computePeriod == 0) {
        Compute(it);
    }
}
void Observable::SetComputePeriod(int iterations) {
    _computePeriod = iterations;
}
// void Observable::SetFilename(const std::string& filename) {
//     _output_filename = filename;
// }

int Observable::MemoryAllocated() const {
    return 0;
}