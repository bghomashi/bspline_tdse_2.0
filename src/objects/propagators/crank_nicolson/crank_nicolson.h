#pragma once

#include "common/objects/propagator.h"
#include "common/maths/math_common.h"
#include "common/maths/gmres_solver.h"
#include "common/maths/vector.h"
#include "common/maths/matrix.h"
#include <complex>
#include <map>

class CrankNicolson : public tdse::Propagator::Register<CrankNicolson> {
    maths::GMRESSolver _solver;

    maths::Vector _psi_temp;
    maths::Matrix _U0p, _U0m, _HI[maths::DimIndex::NUM];
    maths::Matrix _Up, _Um;
public:
    CrankNicolson();

    void Initialize();
    bool DoStep(int it);
    void Finish();

    void FillFieldFree(maths::Matrix& m);
    void FillOverlap(maths::Matrix& m);
    void FillInteractionX(maths::Matrix& m);
    void FillInteractionY(maths::Matrix& m);
    void FillInteractionZ(maths::Matrix& m);
    void FillU0(maths::Matrix& m);

    void DoCheckpoint();
    void DoObservables();

    static std::string GetName();
    static tdse::Propagator::Ptr_t Create(const nlohmann::json& propagator);
    static bool Validate(const nlohmann::json& propagator);
};
