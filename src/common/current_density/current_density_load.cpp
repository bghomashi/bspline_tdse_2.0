#include "common/current_density/current_density.h"
#include "common/system_state/system_state.h"
#include "common/bspline/bspline.h"
#include "common/utility/to_lower.h"
#include "common/maths/matrix.h"

void CurrentDensity::AddPulse(Pulse::Ptr_t pul) {
    _pulses.push_back(pul);
}
void CurrentDensity::Load(const nlohmann::json& input) {
    auto& eigen_state = input["eigen_state"];
    
    // --------- basis -------------
    auto& basis = input["basis"];
    SystemState::Load(input);

    if (input.contains("lasers") && input["lasers"].is_array()) {
        for (auto& pulse : input["lasers"])
            AddPulse(Pulse::Create(pulse["envelope"], pulse));
    }
    _time_step = input["propagator"]["time_step"].get<double>();
}
