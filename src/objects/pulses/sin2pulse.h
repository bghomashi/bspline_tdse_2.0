#ifndef __SIN2PULSE_H__
#define __SIN2PULSE_H__

#include "common/objects/pulse.h"

struct Sin2Pulse : public Pulse::Register<Sin2Pulse> {
    Vec3 A(double t) const;
    Vec3 E(double t) const;

    // static Ptr_t Create(Envelope env, double delay_cycles, double cep, double intensity, double frequency, int numCycles, double ellipticity, const Vec3& polarization, const Vec3& poynting_vector);
    static Ptr_t Create(const nlohmann::json& pusle);
    static bool Validate(const nlohmann::json& pusle);
    static std::string GetName();
};


#endif