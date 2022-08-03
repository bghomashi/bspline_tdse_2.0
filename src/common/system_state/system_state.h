#ifndef __SYSTEM_STATE_H__
#define __SYSTEM_STATE_H__

#include "common/utility/json.hpp"
#include "common/maths/eigen_solver.h"

class SystemState {
    int _lmax_basis, _mmax_basis;
    int _eigenStateBoundNmax, _eigenStateContinuumNmax, _eigenStateBoundLmax, _eigenStateContinuumLmax;
    std::string _eigenStateFilename;
    maths::EigenProblemType _problemType;
public:
    static bool Validate(const nlohmann::json& input);
    static void Load(const nlohmann::json& input);

    static int GetBasisLmax();
    static int GetBasisMmax();
    static int GetEigenStateContinuumLmax();
    static int GetEigenStateBoundLmax();
    static int GetEigenStateContinuumNmax();
    static int GetEigenStateBoundNmax();
    static std::string GetEigenStateFilename();
    static maths::EigenProblemType GetProblemType();
};

#endif