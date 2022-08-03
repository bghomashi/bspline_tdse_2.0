#ifndef __INITIAL_STATE_H__
#define __INITIAL_STATE_H__

#include "common/maths/vector.h"
#include <string>
#include <vector>
#include <memory>
#include "common/utility/json.hpp"


namespace tdse {
    // an object to hold the initial state info
    struct StateDescriptor {
        int n, l, m;
        double phase, amplitude;
    };

    class Simulation;

    class InitialState {
    public:
        typedef std::shared_ptr<InitialState> Ptr_t;
        
        std::vector<StateDescriptor> _initialState;
        std::string _eigenStateFilename;
        int _eigenStateNMax;
        

        void SetFilename(const std::string& filename);
        void BuildInitialState(maths::Vector psi);
        void AddInitialState(int n, int l, int m, double phase,  double amplitude);

        static bool Validate(const nlohmann::json& input);
        void Load(const nlohmann::json& input);
    };
}

#endif