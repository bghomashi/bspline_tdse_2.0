#include "common/objects/pulse.h"
#include "common/utility/logger.h"
#include "common/maths/math_common.h"
#include <sstream>

Pulse::Pulse() : 
    period(0), delay(0), numCycles(0), E0(0), 
    duration(0), frequency(0), intensity(0), ellipticity(0), cep(0),
    polarization_vector(), poynting_vector(), minor_polarization_vector()
{}


