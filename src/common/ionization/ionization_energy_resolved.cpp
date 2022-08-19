#include "common/ionization/ionization.h"
#include "common/maths/math_common.h"
#include "common/utility/spherical_harmonics.h"

void Ionization::EnergyResolved() {
    // HARD CODING NUMBERS FOR TEST
    // energy grid
    std::vector<std::pair<double, double>>  ionization_out(1000);
    double Emin = 0.004, Emax = 1.0, d = (Emax - Emin) / (ionization_out.size() - 1);
    for (int i = 0; i < ionization_out.size(); i++) {
        auto& io = ionization_out[i];

        io.first = Emin + d*i;
        io.second = 0;
    }

}