#include "objects/pulses/sin2pulse.h"
#include "common/maths/math_common.h"
#include "common/utility/logger.h"

Pulse::Ptr_t Sin2Pulse::Create(const nlohmann::json& input) {
    Sin2Pulse* pulse = new Sin2Pulse;

    pulse->numCycles = input["num_cycles"];
    pulse->delay = input["cycles_delay"];
    pulse->intensity = input["intensity"];
    pulse->cep = input["cep"];
    pulse->ellipticity = 0.;
    Vec3 pol_vector, poy_vector;
    double norm = 0;

    pol_vector.x = input["polarization_vector"][0];
    pol_vector.y = input["polarization_vector"][1];
    pol_vector.z = input["polarization_vector"][2];
    pol_vector = normal(pol_vector);
    
    poy_vector.x = input["poynting_vector"][0];
    poy_vector.y = input["poynting_vector"][1];
    poy_vector.z = input["poynting_vector"][2];
    poy_vector = normal(poy_vector);

    // specify wavelength (nm) or energy (au)
    if (input.contains("ellipticity"))
        pulse->ellipticity = input["ellipticity"].get<double>();
    if (input.contains("wavelength"))
        pulse->frequency = maths::LnmToEnergy/input["wavelength"].get<double>();   // not sure why "get..." is needed here
    if (input.contains("energy"))
        pulse->frequency = input["energy"];


    // if (input["envelope"] == "sin2")
    
    pulse->E0 = std::sqrt(pulse->intensity / 3.51e16); 
    pulse->period = 2.*maths::Pi/pulse->frequency;
    pulse->duration = pulse->numCycles*pulse->period;
    pulse->polarization_vector = normal(pol_vector)/std::sqrt(1. + pulse->ellipticity*pulse->ellipticity);
    pulse->poynting_vector = normal(poy_vector);
    pulse->minor_polarization_vector = normal(cross(pol_vector, poy_vector))*(pulse->ellipticity/std::sqrt(1. + pulse->ellipticity*pulse->ellipticity));

    return Pulse::Ptr_t(pulse);
}

bool Sin2Pulse::Validate(const nlohmann::json& pulse) {
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
std::string Sin2Pulse::GetName() {
    return "sin2";
}



Vec3 Sin2Pulse::A(double t) const {     // A(t)
    if (t < delay) return Vec3{0};

    double T = t-delay+cep;
    double Env = -E0*sin(maths::Pi*T/duration)*sin(maths::Pi*T/duration)/frequency;
    Vec3 p = sin(frequency*T)*polarization_vector - cos(frequency*T)*minor_polarization_vector;
    return Env*p;
}

Vec3 Sin2Pulse::E(double t) const {     // E=-dA/dt
    if (t < delay) return Vec3{0};

    double T = t-delay+cep;
    // product rule
    double Env1 = E0*sin(maths::Pi*T/duration)*sin(maths::Pi*T/duration);                             // dont diff. env
    double Env2 = 2.*maths::Pi*E0*sin(maths::Pi*T/duration)*cos(maths::Pi*T/duration)/frequency/duration;    // do diff. env

    Vec3 p1 = cos(frequency*T)*polarization_vector + sin(frequency*T)*minor_polarization_vector;       // do diff. carrier wave
    Vec3 p2 = sin(frequency*T)*polarization_vector - cos(frequency*T)*minor_polarization_vector;                           // dont

    return Env1*p1 + Env2*p2;

}