
#include "common/tdse/simulation.h"
#include "common/utility/logger.h"

using namespace tdse;

static tdse::Simulation s_simulation;

Simulation::Simulation() : 
    _do_propagate(true), _restarting(false), 
    _checkpoints(0), _laserAxial(true), 
    _potentialSymmetry(Symmetry::Central) {
    _pol[maths::X] = _pol[maths::Y] = _pol[maths::Z] = false;
}

int Simulation::GetDOF() {
    return s_simulation._dof;
}
double Simulation::GetTmax() {
    return s_simulation._tmax;
}
double Simulation::GetTimestep() {
    return s_simulation._propagator->GetTimeStep(); //_dt;
}
double Simulation::GetTime() {
    return s_simulation._propagator->GetLastTime();
}
int Simulation::GetNumTimeSteps() {
    return s_simulation._NT;
}
const bool* Simulation::GetPolarization() {
    return s_simulation._pol;
}
const std::vector<double>& Simulation::GetField(int dim_index) {
    return s_simulation._field[dim_index];
}
const std::vector<double>* Simulation::GetFields() {
    return s_simulation._field;
}
const std::vector<Pulse::Ptr_t>& Simulation::GetPulses() {
    return s_simulation._pulses;
}
const std::vector<int>& Simulation::GetMs() {
    return s_simulation._Ms;
}
const std::vector<int>& Simulation::GetMRows() {
    return s_simulation._mRows;
}
const std::vector<Potential::Ptr_t>& Simulation::GetPotentials() {
    return s_simulation._potentials;
}
maths::Vector& Simulation::GetPsi() {
    return s_simulation._psi;
}
std::string Simulation::GetInitialStateFile() {
    return s_simulation._initialState._eigenStateFilename;
}
int Simulation::GetInitialStateNmax() {
    return s_simulation._initialState._eigenStateNMax;
}
int Simulation::GetPotentialSymmetry() {
    return s_simulation._potentialSymmetry;
}
bool Simulation::GetLaserSymmetry() {
    return s_simulation._laserAxial;
}




void Simulation::AddPotential(Potential::Ptr_t pot) {
    _potentials.push_back(pot);
}
void Simulation::AddPulse(Pulse::Ptr_t pul) {
    _pulses.push_back(pul);
}
void Simulation::AddObservable(Observable::Ptr_t obs) {
    _observables.push_back(obs);   
}


void Simulation::Load(const nlohmann::json& input) {
    s_simulation._Load(input);
}
bool Simulation::Validate(const nlohmann::json& input) {
    return s_simulation._Validate(input);
}

void Simulation::Execute() {
    s_simulation._Execute();
}