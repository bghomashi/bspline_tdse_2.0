
// utility
// read in input
// read in angularly resolved ionization output
// read in eigenstates
// bspline basis
// A) output energy resolved results 
// B) output density
#include <iostream>
#include <fstream>
#include <limits>
#include <unordered_map>

#include "common/file_io/io_factory.h"
#include "common/system_state/system_state.h"
#include "common/bspline/bspline.h"
#include "common/maths/math_factory.h"
#include "common/utility/logger.h"
#include "common/utility/profiler.h"
#include "common/utility/to_lower.h"
#include "common/utility/json.hpp"
#include "petsc/maths/petsc_math_factory.h"
#include "petsc/utility/petsc_logger.h"
#include "petsc/utility/petsc_profiler.h"
#include "petsc/file_io/petsc_io_factory.h"


bool ValidateAndRunIonization(int argc, char **args, const std::string& filename);

int main(int argc, char **args) {
    if (!ValidateAndRunIonization(argc, args,"input.json"))
        return -1;
    return 0;
}

class Ionization {
    bool _ecs_on;
    std::string _ionizationFilename ;
public:
    bool Validate(const nlohmann::json& input) {
        LOG_INFO("Validating input file.");

        if (!SystemState::Validate(input))
            return false;
        
        auto& eigen_state = input["eigen_state"];
        
        // --------- ecs
        if (eigen_state.contains("ecs_on") && !eigen_state["ecs_on"].is_boolean()) {
            LOG_CRITICAL("Optional entry \"ecs_on\" in eigen_state must be a boolean.");
            return false; 
        }
        // ---------------- basis ----------------
        if (!bspline::Basis::Validate(input["basis"])) {
            LOG_CRITICAL("Failed to contruct \"basis\".");
            return false;
        }

        // --------------- ionization observable ---------
        if (!input.contains("observables")) {
            LOG_CRITICAL("Input file requires \"observables\" object.");
            return false;
        }
        if (!input["observables"].contains("ionization")) {
            LOG_CRITICAL("\"Observables\" requires object entry: ionization.");
            return false;
        }
        if (!input["observables"]["ionization"].contains("filename") || !input["observables"]["ionization"]["filename"].is_string()) {
            LOG_CRITICAL("\"ionization\" requires string entry: filename.");
            return false;
        }
        
        return true;
    }

    void Load(const nlohmann::json& input) {
        auto& eigen_state = input["eigen_state"];
        // --------- ecs -----------------
        if (eigen_state.contains("ecs_on"))
            _ecs_on = eigen_state["ecs_on"].get<bool>();

        // --------- basis -------------
        auto& basis = input["basis"];
        SystemState::Load(input);
        bspline::Basis::Load(input["basis"], !_ecs_on);

        // --------- ionization filename -------
        _ionizationFilename = input["observables"]["ionization"]["filename"].get<std::string>();
    }

    void Execute() {
        std::unordered_map<int,                 // l
            std::unordered_map<int,             // m
                std::vector<std::pair<double,   // energy
                    double>>>>                  // population
                 populations;

        // HARD CODING NUMBERS FOR TEST
        std::vector<std::pair<double, double>>  ionization_out(1000);
        double Emin = 0.004, Emax = 0.15, d = (Emax - Emin) / (ionization_out.size() - 1);
        for (int i = 0; i < ionization_out.size(); i++) {
            auto& io = ionization_out[i];

            io.first = Emin + d*i;
            io.second = 0;
        }

        double energy, amp_real, amp_imag, pop;
        int l, m;
        std::ifstream ionization_in(_ionizationFilename);
        ionization_in.ignore(std::numeric_limits<std::streamsize>::max(), ionization_in.widen('\n'));
        
        while (!ionization_in.eof()) {
            ionization_in >> energy >> l >> m >> pop >> amp_real >> amp_imag;
            populations[l][m].push_back({energy, pop});
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
        for (auto& lPop : populations) {
            int l = lPop.first;
            for (auto& mPop : lPop.second) {
                int m = mPop.first;

                auto& energyVector = mPop.second;
                for (int n = 0; n < energyVector.size()-1; n++) {
                    double E1 = energyVector[n].first, E2 = energyVector[n+1].first;
                    double dE = E2 - E1;
                    double pop = energyVector[n].second/dE;

                    auto interval = findInterval(E1, E2);
                    for (int i = interval.first; i <= interval.second; i++)
                        ionization_out[i].second += pop;
                }
            }
        }

        std::ofstream outFile("energy_resolved.txt");
        outFile << std::setprecision(8) << std::scientific;
        
        for (auto& EPop : ionization_out) {
            double E = EPop.first;
            double pop = EPop.second;

            outFile << E << "\t" << pop << std::endl;
        }
    }
};


bool ValidateAndRunIonization(int argc, char **args, const std::string& filename) {
     std::ifstream i(filename);
    nlohmann::json input;
    i >> input;

    // log filename is optional. If not found output to stdout
    if (input.contains("log_filename") && input["log_filename"].is_string())
        Log::SetLoggerFile(input["log_filename"]);

    // first we immediately check the math library
    // - if this is library has a special 'logger' we start it right now
    // if (!ValidateMathLibrary(input))
    //     return false;

    if (ToLower(input["math_library"]) == "petsc") {
        maths::Factory::SetInstance(new PetscMathFactory());
        io::Factory::SetInstance(new PetscIOFactory());
        Log::SetLogger(new PetscLogger());
        Profile::SetProfiler(new PetscProfiler());
    } else if (ToLower(input["math_library"]) == "thread_pool") {
        std::cout << "thread_pool is not yet supported" << std::endl;
        return false;
    }

    if (!maths::Factory::Startup(argc, args))
        return false;
    if (!io::Factory::Startup())
        return false;

    Ionization ionCalc;
    if (!ionCalc.Validate(input))
        return false;

    ionCalc.Load(input);
    ionCalc.Execute();

    
    LOG_INFO("Shutting down.\n------------------------------------------------\n\n");
    Profile::PrintTo("profile.txt");
    io::Factory::Shutdown();
    maths::Factory::Shutdown();

    return true;
}

