#include "common/populations/populations.h"

#include "common/system_state/system_state.h"
#include "common/bspline/bspline.h"
#include "common/utility/to_lower.h"
#include "common/maths/matrix.h"

void Populations::Load(const nlohmann::json& input) {
    auto& eigen_state = input["eigen_state"];
    
    // --------- ecs -----------------
    if (eigen_state.contains("ecs_on"))
        _ecs_on = eigen_state["ecs_on"].get<bool>();

    // --------- basis -------------
    auto& basis = input["basis"];
    SystemState::Load(input);
    bspline::Basis::Load(input["basis"], !_ecs_on);

    _eigenstate_filename = eigen_state["filename"];

    // --------- populations filename -------
    _populations_filename = input["observables"]["populations"]["filename"].get<std::string>();
    
}
