#include "objects/pulses/trapezoidal_pulse.h"
#include "common/maths/math_common.h"
#include "common/utility/logger.h"

Pulse::Ptr_t TrapezoidalPulse::Create(const nlohmann::json& input) {
    TrapezoidalPulse* pulse = new TrapezoidalPulse;

    pulse->numCycles = input["num_cycles"].get<int>();
    pulse->delay = input["cycles_delay"].get<double>();
    pulse->intensity = input["intensity"].get<double>();
    pulse->cep = input["cep"].get<double>()*2.*maths::Pi;
    pulse->ellipticity = 0.;
    Vec3 pol_vector, poy_vector;
    double norm = 0;

    pol_vector.x = input["polarization_vector"][0].get<double>();
    pol_vector.y = input["polarization_vector"][1].get<double>();
    pol_vector.z = input["polarization_vector"][2].get<double>();
    pol_vector = normal(pol_vector);
    
    poy_vector.x = input["poynting_vector"][0].get<double>();
    poy_vector.y = input["poynting_vector"][1].get<double>();
    poy_vector.z = input["poynting_vector"][2].get<double>();
    poy_vector = normal(poy_vector);

    pulse->cycles_up = input["cycles_up"].get<double>();
    pulse->cycles_down = input["cycles_down"].get<double>();
    pulse->cycles_plateau = pulse->numCycles - pulse->cycles_down - pulse->cycles_up;

    // specify wavelength (nm) or energy (au)
    if (input.contains("ellipticity"))
        pulse->ellipticity = input["ellipticity"].get<double>();
    if (input.contains("wavelength"))
        pulse->frequency = maths::LnmToEnergy/input["wavelength"].get<double>();   // not sure why "get..." is needed here
    if (input.contains("energy"))
        pulse->frequency = input["energy"];


    pulse->E0 = std::sqrt(pulse->intensity / 3.51e16); 
    pulse->period = 2.*maths::Pi/pulse->frequency;
    pulse->duration = pulse->numCycles*pulse->period;
    pulse->polarization_vector = normal(pol_vector)/std::sqrt(1. + pulse->ellipticity*pulse->ellipticity);
    pulse->poynting_vector = normal(poy_vector);
    pulse->minor_polarization_vector = normal(cross(pol_vector, poy_vector))*(pulse->ellipticity/std::sqrt(1. + pulse->ellipticity*pulse->ellipticity));

    return Pulse::Ptr_t(pulse);
}

bool TrapezoidalPulse::Validate(const nlohmann::json& pulse) {
    if (!(pulse.contains("cycles_delay") && pulse["cycles_delay"].is_number())) {
        LOG_CRITICAL("each \"laser\" must contain number entry: cycles_delay");
        return false;
    }
    if (!(pulse.contains("envelope") && pulse["envelope"].is_string())) {
        LOG_CRITICAL("each \"laser\" must contain string entry: envelope");
        return false;
    }
    if (!(pulse.contains("intensity") && pulse["intensity"].is_number())) {
        LOG_CRITICAL("each \"laser\" must contain number entry: intensity");
        return false;
    }
    if (!(pulse.contains("num_cycles") && pulse["num_cycles"].is_number())) {
        LOG_CRITICAL("each \"laser\" must contain number entry: num_cycles");
        return false;
    }
    if (!(pulse.contains("wavelength") && pulse["wavelength"].is_number()) && 
        !(pulse.contains("energy") && pulse["energy"].is_number())) {
        LOG_CRITICAL("each \"laser\" must contain number entry: either wavelength or energy.");
        return false;
    }
    if (!(pulse.contains("cep") && pulse["cep"].is_number())) {
        LOG_CRITICAL("each \"laser\" must contain number entry: cep");
        return false;
    }
    if (!(pulse.contains("cycles_up") && pulse["cycles_up"].is_number())) {
        LOG_CRITICAL("\"trapezoidal\" laser must contain number entry: cycles_up");
        return false;
    }
    if (!(pulse.contains("cycles_down") && pulse["cycles_down"].is_number())) {
        LOG_CRITICAL("\"trapezoidal\" laser must contain number entry: cycles_down");
        return false;
    }
    if (pulse.contains("ellipticity") && !pulse["ellipticity"].is_number()) {
        LOG_CRITICAL("optional parameter \"ellipticity\" must be a number.");
        return false;
    }
    if (!(pulse.contains("polarization_vector") && pulse["polarization_vector"].is_array())) {
        LOG_CRITICAL("each \"laser\" must contain 3-vector entry: polarization_vector.");
        return false;
    }
    if (!(pulse.contains("poynting_vector") && pulse["poynting_vector"].is_array())) {
        LOG_CRITICAL("each \"laser\" must contain 3-vector entry: poynting_vector.");
        return false;
    }
    return true;
}
std::string TrapezoidalPulse::GetName() {
    return "trapezoidal";
}


double envelope(double t, double cycles_up, double cycles_down, double cycles_plateau, double frequency) {
    double t1 = cycles_up*2.*maths::Pi/frequency;
    double t2 = cycles_plateau*2.*maths::Pi/frequency + t1;
    double t3 = cycles_down*2.*maths::Pi/frequency + t2;
    if (t <= 0)
        return 0.0;
    else if (t <= t1)
        return t/t1;
    else if (t <= t2)
        return 1.0;
    else if (t <= t3)
        return (t3-t) / (t3-t2);
    return 0.0;
}
double d_envelope(double t, double cycles_up, double cycles_down, double cycles_plateau, double frequency) {
    double t1 = cycles_up*2.*maths::Pi/frequency;
    double t2 = cycles_plateau*2.*maths::Pi/frequency + t1;
    double t3 = cycles_down*2.*maths::Pi/frequency + t2;
    if (t <= 0)
        return 0.0;
    else if (t <= t1)
        return 1./t1;
    else if (t <= t2)
        return 0.0;
    else if (t <= t3)
        return -1. / (t3-t2);
    return 0.0;
};

Vec3 TrapezoidalPulse::A(double t) const {     // A(t)
    if (t < delay) return Vec3{0};

    double T = t-delay;
    double Env = -E0*envelope(T, cycles_up, cycles_down, cycles_plateau, frequency);
    Vec3 p = sin(frequency*T+cep)*polarization_vector - cos(frequency*T+cep)*minor_polarization_vector;
    return Env*p;
}

Vec3 TrapezoidalPulse::E(double t) const {     // E=-dA/dt
    if (t < delay) return Vec3{0};

    double T = t-delay;
    // product rule
    double Env1 = E0*envelope(T, cycles_up, cycles_down, cycles_plateau, frequency);                // dont diff. env
    double Env2 = E0*d_envelope(T, cycles_up, cycles_down, cycles_plateau, frequency)/frequency;    // do diff. env

    Vec3 p1 = cos(frequency*T+cep)*polarization_vector + sin(frequency*T+cep)*minor_polarization_vector;       // do diff. carrier wave
    Vec3 p2 = sin(frequency*T+cep)*polarization_vector - cos(frequency*T+cep)*minor_polarization_vector;                           // dont

    return Env1*p1 + Env2*p2;

}