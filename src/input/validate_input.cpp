#include "input/validate_input.h"
#include <iostream>
#include <fstream>
#include "common/maths/math_factory.h"
#include "common/utility/json.hpp"
#include "common/utility/logger.h"
#include "common/file_io/io_factory.h"
#include "petsc/maths/petsc_math_factory.h"
#include "petsc/file_io/petsc_io_factory.h"
#include "petsc/utility/petsc_logger.h"
#include "petsc/utility/petsc_profiler.h"
#include "common/utility/to_lower.h"
#include "common/tdse/simulation.h"
#include "common/tise/tise.h"

bool ValidateAndRunTDSE(int argc, char **args, const std::string& filename) {
    std::ifstream i(filename);
    nlohmann::json input;
    i >> input;

    // log filename is optional. If not found output to stdout
    if (input.contains("log_filename") && input["log_filename"].is_string())
        Log::SetLoggerFile(input["log_filename"]);

    // first we immediately check the math library
    // - if this is library has a special 'logger' we start it right now
    // if (!ValidateMathLibrary(input))
    //     return false;

    if (ToLower(input["math_library"]) == "petsc") {
        maths::Factory::SetInstance(new PetscMathFactory());
        io::Factory::SetInstance(new PetscIOFactory());
        Log::SetLogger(new PetscLogger());
        Profile::SetProfiler(new PetscProfiler());
    } else if (ToLower(input["math_library"]) == "thread_pool") {
        std::cout << "thread_pool is not yet supported" << std::endl;
        return false;
    }

    if (!maths::Factory::Startup(argc, args))
        return false;
    if (!io::Factory::Startup())
        return false;
    if (!tdse::Simulation::Validate(input))
        return false;


    tdse::Simulation::Load(input);
    tdse::Simulation::Execute();

    LOG_INFO("Shutting down.\n------------------------------------------------\n\n");
    Profile::PrintTo("profile.txt");
    io::Factory::Shutdown();
    maths::Factory::Shutdown();

    return true;
}

bool ValidateRunTISE(int argc, char **args, const std::string& filename) {
    std::ifstream i(filename);
    nlohmann::json input;
    i >> input;

    // log filename is optional. If not found output to stdout
    if (input.contains("log_filename") && input["log_filename"].is_string())
        Log::SetLoggerFile(input["log_filename"]);

    // first we immediately check the math library
    // - if this is library has a special 'logger' we start it right now
    // if (!ValidateMathLibrary(input))
    //     return false;

    if (ToLower(input["math_library"]) == "petsc") {
        maths::Factory::SetInstance(new PetscMathFactory());
        io::Factory::SetInstance(new PetscIOFactory());
        Log::SetLogger(new PetscLogger());
        Profile::SetProfiler(new PetscProfiler());
    } else if (ToLower(input["math_library"]) == "thread_pool") {
        std::cout << "thread_pool is not yet supported" << std::endl;
        return false;
    }

    if (!maths::Factory::Startup(argc, args))
        return false;
    if (!io::Factory::Startup())
        return false;
    if (!tise::TISE::Validate(input))
        return false;

    tise::TISE::Load(input);
    tise::TISE::Execute();

    LOG_INFO("Shutting down.\n------------------------------------------------\n\n");
    Profile::PrintTo("profile.txt");
    io::Factory::Shutdown();
    maths::Factory::Shutdown();

    return true;
}