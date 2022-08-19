#ifndef __TRAPEZOIDAL_PULSE_H__
#define __TRAPEZOIDAL_PULSE_H__

#include "common/objects/pulse.h"

struct TrapezoidalPulse : public Pulse::Register<TrapezoidalPulse> {
    double cycles_up, cycles_down, cycles_plateau;

    Vec3 A(double t) const;
    Vec3 E(double t) const;

    // static Ptr_t Create(Envelope env, double delay_cycles, double cep, double intensity, double frequency, int numCycles, double ellipticity, const Vec3& polarization, const Vec3& poynting_vector);
    static Ptr_t Create(const nlohmann::json& pusle);
    static bool Validate(const nlohmann::json& pusle);
    static std::string GetName();
};


#endif