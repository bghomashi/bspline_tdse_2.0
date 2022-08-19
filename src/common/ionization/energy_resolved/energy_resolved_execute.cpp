#include "common/ionization/ionization.h"
#include "common/ionization/energy_resolved/energy_resolved.h"
#include <iostream>
#include <fstream>
#include <iomanip>

void EnergyResolved::Execute() {
    auto& populations = _parent->populations;

    // energy grid
    int n = (_emax - _emin) / (_estep + 1);
    std::vector<std::pair<double, double>>  ionization_out(n);
    for (int i = 0; i < ionization_out.size(); i++) {
        auto& io = ionization_out[i];

        io.first = _emin + _estep*i;
        io.second = 0;
    }


    auto findInterval = [&](double E1, double E2) {
        int start = 0, end = 0;
        // find first energy >= E1
        for (start = 0; start < ionization_out.size(); start++) {
            if (ionization_out[start].first >= E1)
                break;
        }
        for (end = ionization_out.size()-1; end >= 0 ; end--) {
            if (ionization_out[end].first <= E2)
                break;
        }
        return std::pair<int, int>(start, end);
    };

    // we want to sum over all the L's and M's
    // for (auto& lPop : populations)
    for (int l = 0; l < 3; l++)
    {
        // int l = lPop.first;
        auto& ms = populations[l];
        for (auto& mPop : ms) {
            int m = mPop.first;

            auto& energyVector = mPop.second;
            for (int n = 0; n < energyVector.size()-1; n++) {
                double E1 = energyVector[n].energy, E2 = energyVector[n+1].energy;
                double dE = E2 - E1;
                double pop = energyVector[n].population/dE;

                auto interval = findInterval(E1, E2);
                for (int i = interval.first; i <= interval.second; i++)
                    ionization_out[i].second += pop;
            }
        }
    }

    std::ofstream outfile(_output_filename);
    outfile << std::setprecision(8) << std::scientific;
    
    for (auto& EPop : ionization_out) {
        double E = EPop.first;
        double pop = EPop.second;

        outfile << E << "\t" << pop << std::endl;
    }

}