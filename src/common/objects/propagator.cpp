#include "common/objects/propagator.h"
#include "common/tdse/simulation.h"

using namespace tdse;

Propagator::Propagator() : _t(0), _dt(0) {}

double Propagator::GetLastTime() const {
    return _t;
}
double Propagator::GetTimeStep() const {
    return _dt;
}