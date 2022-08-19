#include "common/ionization/ionization.h"

#include "common/system_state/system_state.h"
#include "common/bspline/bspline.h"
#include "common/utility/to_lower.h"

void Ionization::Load(const nlohmann::json& input) {
    auto& eigen_state = input["eigen_state"];
    auto& ionization = input["ionization"];
    // --------- ecs -----------------
    if (eigen_state.contains("ecs_on"))
        _ecs_on = eigen_state["ecs_on"].get<bool>();

    // --------- basis -------------
    auto& basis = input["basis"];
    SystemState::Load(input);
    bspline::Basis::Load(input["basis"], !_ecs_on);

    _eigenstate_filename = eigen_state["filename"];
    continuum_nmax = eigen_state["continuum"]["nmax"];

    // --------- ionization filename -------
    _ionization_filename = input["observables"]["ionization"]["filename"].get<std::string>();
    
    // --------- ionization ------------
    for (auto& iot_pair : ionization.items()) {
        _ionization_tasks.push_back(IonizationTask::Create(iot_pair.key(), iot_pair.value()));
        _ionization_tasks.back()->_parent = this;
    }

}
