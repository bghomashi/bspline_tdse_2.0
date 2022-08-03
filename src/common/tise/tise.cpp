#include "common/tise/tise.h"

using namespace tise;


static tise::TISE s_tise;

TISE::TISE() : _tol(1e-10), _ecs_on(false) {}


void TISE::AddPotential(tdse::Potential::Ptr_t p) {
    _potentials.push_back(p);
}
void TISE::_Execute() {
    Solve();

    _eigensolver = nullptr;
    _tise_out = nullptr;
    for (auto& p : _potentials)
        p = nullptr;
}


void TISE::Execute() {
    s_tise._Execute();
}
bool TISE::Validate(const nlohmann::json& input) {
    return s_tise._Validate(input);
}
void TISE::Load(const nlohmann::json& input) {
    s_tise._Load(input);
}