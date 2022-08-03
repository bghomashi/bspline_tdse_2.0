#ifndef __PARAMETERS_H__
#define __PARAMETERS_H__

#include "common/utility/objects/json.hpp"

namespace tdse {
    struct Parameters {
        // radial basis parameters
        double rmin, rmax;
        int order, nodes;
        double ecs_r0, ecs_theta;
        bool skip_first, skip_last;

        // angular basis parameters
        int lmax, mmax;

        // temporal parameters
        double dt;

        // checkpoints
        int checkpointIterationPeriod;
        bool restartingFromCheckpoint, doPropagate;

        static bool ValidateAndStoreBasis(const nlohmann::json& input);
        static bool ValidateAndStoreSimulation(const nlohmann::json& input);
    };
}

#endif