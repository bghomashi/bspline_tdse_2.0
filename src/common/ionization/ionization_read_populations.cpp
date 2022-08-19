#include "common/ionization/ionization.h"
#include <iostream>
#include <fstream>
#include <climits>

void Ionization::ReadPopulations() {
    double energy, amp_real, amp_imag, pop, phase;
    int l, m;
    std::ifstream ionization_in(_ionization_filename);
    ionization_in.ignore(std::numeric_limits<std::streamsize>::max(), ionization_in.widen('\n'));
    
    while (!ionization_in.eof()) {
        ionization_in >> energy >> l >> m >> pop >> amp_real >> amp_imag >> phase;
        if (ionization_in.fail()) break;
        populations[l][m].push_back({energy, pop, phase, {amp_real, amp_imag}});
    }


    lmax = INT_MIN, mmax = INT_MIN;

    for (auto& l : populations) {
        if (l.first > lmax)               // compare l values
            lmax = l.first;

        for (auto& m : l.second) {
            if (m.first > mmax)
                mmax = m.first;
        }
    }

    // PROBABLY NOT GOOD:
    // this establishes a sign convention
    // maths::complex ddE0, ddE1, ddE2;
    // for (auto& l : populations) {
    //     for (auto& m : l.second) {
    //         auto& amps = m.second;
    //         for (int i = 2; i < amps.size()-1; i++) {
    //             // backward difference second derivative
    //             ddE0 = amps[i-2].amplitude - 2.*amps[i-1].amplitude + amps[i].amplitude /
    //                     ((amps[i].energy - amps[i-1].energy)*(amps[i-1].energy - amps[i-2].energy));

    //             // center difference 1
    //             ddE1 = (amps[i+1].amplitude - 2.*amps[i].amplitude + amps[i-1].amplitude) / 
    //                     ((amps[i+1].energy - amps[i].energy)*(amps[i].energy - amps[i-1].energy));

    //             // center difference 2
    //             ddE2 = (-amps[i+1].amplitude - 2.*amps[i].amplitude + amps[i-1].amplitude) / 
    //                     ((amps[i+1].energy - amps[i].energy)*(amps[i].energy - amps[i-1].energy));

    //             if (std::abs(ddE1 - ddE0) > std::abs(ddE2 - ddE0)) {
    //                 amps[i+1].amplitude *= -1;
    //             }

    //         }
    //     }
    // }




}