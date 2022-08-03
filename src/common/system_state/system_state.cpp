#include "common/system_state/system_state.h"
#include "common/utility/logger.h"
#include "common/utility/file_exists.h"
#include "common/utility/to_lower.h"

static SystemState s_systemState;



bool SystemState::Validate(const nlohmann::json& input) {
    
    if (!(input.contains("eigen_state") && input["eigen_state"].is_object())) {
        LOG_CRITICAL("input file must contain object entry: eigen_state");
        return false;
    }
    // -------------------- validate eigenstate input -------------------
    auto& eigenState = input["eigen_state"];
    if (!(eigenState.contains("tol") && eigenState["tol"].is_number())) {
        LOG_CRITICAL("\"eigen_state\" must contain number entry: tol");
        return false;
    }
    if (eigenState.contains("problem_type") && !eigenState["problem_type"].is_string()) {
        LOG_CRITICAL("Optional entry \"problem_type\" in eigen_state must be a string.");
        return false;
    }
    if (eigenState.contains("problem_type")) {
        if (ToLower(eigenState["problem_type"].get<std::string>()) != "ghep" && 
            ToLower(eigenState["problem_type"].get<std::string>()) != "gnhep" &&
            ToLower(eigenState["problem_type"].get<std::string>()) != "pgnhep") {
            LOG_CRITICAL("\"problem_type\" must be one of: \"GHEP\", \"GNHEP\", \"PGNHEP\".");
            return false;
        }
    }

    // -------------------- filename  ---------------
    if (!(eigenState.contains("filename") && eigenState["filename"].is_string())) {
        LOG_CRITICAL("\"eigen_state\" must contain string entry: filename");
        return false;
    }

    if (!(eigenState.contains("bound") && eigenState["bound"].is_object())) {
        LOG_CRITICAL("\"eigen_state\" must contain object entry: bound");
        return false;
    }
    
    auto& bound = eigenState["bound"];
    // -------------------- bound nmax -----------------
    if (!bound.contains("nmax") || !bound["nmax"].is_number()) {
        LOG_CRITICAL("\"bound\" must contain number entry: nmax.");
        return false; 
    }
    if (bound["nmax"].get<int>() <= 0) {
        LOG_CRITICAL("\"bound\" entry \"nmax\" must be positive.");
        return false; 
    }
    // -------------------- bound lmax -----------------
    if (bound.contains("lmax") && !bound["lmax"].is_number()) {
        LOG_CRITICAL("Optional entry \"lmax\" in bound must be a number.");
        return false; 
    }
    if (bound.contains("lmax") && bound["lmax"].get<int>() < 0) {
        LOG_CRITICAL("bound \"lmax\" must be positive or zero.");
        return false; 
    }
    if (bound.contains("lmax") && bound["lmax"].get<int>() > bound["nmax"].get<int>()) {
        LOG_CRITICAL("bound requires \"lmax\" <= \"nmax\".");
        return false; 
    }

    if (eigenState.contains("continuum")) {
        auto& continuum = eigenState["continuum"];
        // -------------------- continuum nmax -----------------
        if (!continuum.contains("nmax") || !continuum["nmax"].is_number()) {
            LOG_CRITICAL("\"continuum\" must contain number entry: nmax.");
            return false; 
        }
        if (continuum["nmax"].get<int>() <= 0) {
            LOG_CRITICAL("\"continuum\" entry \"nmax\" must be positive.");
            return false; 
        }
        // -------------------- continuum lmax -----------------
            if (!continuum.contains("lmax") || !continuum["lmax"].is_number()) {
            LOG_CRITICAL("\"continuum\" must contain number entry: lmax.");
            return false; 
        }
        if (continuum["lmax"].get<int>() < 0) {
            LOG_CRITICAL("continuum \"lmax\" must be positive or zero.");
            return false; 
        }
    }

    // ----------- has basis -----------
    if (!input.contains("basis") || !input["basis"].is_object()) {
        LOG_CRITICAL("input file must contain object entry: basis");
        return false;
    }

    auto& basis = input["basis"];
    if (!(basis.contains("lmax") && basis["lmax"].is_number())) {
        LOG_CRITICAL("basis must contain number entry: lmax");
        return false;
    }
    if (!(basis.contains("mmax") && basis["mmax"].is_number())) {
        LOG_CRITICAL("basis must contain number entry: mmax");
        return false;
    }
    if (basis["mmax"] > basis["lmax"]) {
        LOG_CRITICAL("mmax must be less than or equal to lmax.");
        return false;
    }
    return true;
}
void SystemState::Load(const nlohmann::json& input) {
    auto& eigenState = input["eigen_state"];
    // ------------ basis ----------------
    s_systemState._lmax_basis = input["basis"]["lmax"].get<int>();
    s_systemState._mmax_basis = input["basis"]["mmax"].get<int>();
    // ------------ defaults --------------- 
    s_systemState._eigenStateBoundLmax = 0;
    s_systemState._eigenStateBoundNmax = 0;
    s_systemState._eigenStateContinuumNmax = 0;
    s_systemState._eigenStateContinuumLmax = 0;
    s_systemState._problemType = maths::EigenProblemType::GNHEP;
    s_systemState._eigenStateFilename = eigenState["filename"].get<std::string>();
    if (eigenState.contains("problem_type")) {
        if (ToLower(eigenState["problem_type"].get<std::string>()) == "ghep")
            s_systemState._problemType = maths::EigenProblemType::GHEP;
        else if (ToLower(eigenState["problem_type"].get<std::string>()) == "gnhep")
            s_systemState._problemType = maths::EigenProblemType::GNHEP;
        else if (ToLower(eigenState["problem_type"].get<std::string>()) == "pgnhep")
            s_systemState._problemType = maths::EigenProblemType::PGNHEP;
    }
    // ------------ bound states ----------- (bound states are required)
    s_systemState._eigenStateBoundNmax = eigenState["bound"]["nmax"].get<int>();
    s_systemState._eigenStateBoundLmax = s_systemState._eigenStateBoundNmax-1;

    if (eigenState["bound"].contains("lmax"))
        s_systemState._eigenStateBoundLmax = eigenState["bound"]["lmax"].get<int>();

    // ------------ continuum states ----------- 
    if (eigenState.contains("continuum")) {
        s_systemState._eigenStateContinuumNmax = eigenState["continuum"]["nmax"].get<int>();
        s_systemState._eigenStateContinuumLmax = eigenState["continuum"]["lmax"].get<int>();
    }
}

int SystemState::GetBasisLmax() {
    return s_systemState._lmax_basis;
}
int SystemState::GetBasisMmax() {
    return s_systemState._mmax_basis;
}
int SystemState::GetEigenStateBoundLmax() {
    return s_systemState._eigenStateBoundLmax;
}
int SystemState::GetEigenStateContinuumLmax() {
    return s_systemState._eigenStateContinuumLmax;
}
int SystemState::GetEigenStateBoundNmax() {
    return s_systemState._eigenStateBoundNmax;
}
int SystemState::GetEigenStateContinuumNmax() {
    return s_systemState._eigenStateContinuumNmax;
}
std::string SystemState::GetEigenStateFilename() {
    return s_systemState._eigenStateFilename;
}
maths::EigenProblemType SystemState::GetProblemType() {
    return s_systemState._problemType;
}