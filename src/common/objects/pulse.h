#ifndef __PULSE_H__
#define __PULSE_H__

#include "common/utility/vec3.h"
#include "common/utility/factory.h"
#include <memory>
#include <string>

struct Pulse : public utl::Factory<Pulse> {
public:
    typedef std::shared_ptr<Pulse> Ptr_t;

    int numCycles;
    double E0, frequency, intensity;
    double period, duration, delay, cep;
    double ellipticity;
    Vec3 polarization_vector;
    Vec3 poynting_vector;
    Vec3 minor_polarization_vector;
    
    Pulse();


    virtual Vec3 A(double t) const = 0;
    virtual Vec3 E(double t) const = 0;


};





// GaussianPulse
// BoxPulse

#endif