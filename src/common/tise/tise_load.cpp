#include "common/tise/tise.h"
#include "common/utility/logger.h"
#include "common/bspline/bspline.h"
#include "common/system_state/system_state.h"

using namespace tise;

void TISE::_Load(const nlohmann::json& input) {
    auto& eigen_state = input["eigen_state"];
    // --------- ecs -----------------
    if (eigen_state.contains("ecs_on"))
        _ecs_on = eigen_state["ecs_on"].get<bool>();
    if (eigen_state.contains("tol"))
        _tol = eigen_state["tol"].get<double>();

    // --------- basis -------------
    auto& basis = input["basis"];
    bspline::Basis::Load(input["basis"], !_ecs_on);
    SystemState::Load(input);

    // --------- solver -------------
    _eigensolver = maths::Factory::CreateEigenSolver();

    // --------- potentials -------------
    if (input.contains("potentials")) {
        for (auto& pot : input["potentials"]) {
            AddPotential(tdse::Potential::Create(pot["type"], pot));
        }
    }
}